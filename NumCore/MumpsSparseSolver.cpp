/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
 the City of New York, and others.
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.*/

#include "stdafx.h"

#include <cassert>
#include <cstring>
#include <vector>

#include "MatrixTools.h"
#include <FECore/CompactUnSymmMatrix.h>
#include <FECore/CompactSymmMatrix.h>
#include <FECore/log.h>

#include "MumpsSparseSolver.h"

#ifdef HAVE_MUMPS
extern "C" {
#include <dmumps_c.h>
}
#ifndef JOB_INIT
#define JOB_INIT -1
#endif
#ifndef JOB_END
#define JOB_END -2
#endif
#endif

class MumpsSparseSolver::Impl
{
public:
    FEModel* fem = nullptr;
    CompactMatrix* A = nullptr;

    int printLevel = 0;
    int matrixMode = MM_AUTO;
    bool isSymmetric = false;
    bool isFactored = false;
    int n = 0;
    int nnz = 0;

#ifdef HAVE_MUMPS
    DMUMPS_STRUC_C id;
    bool initialized = false;

    // MUMPS wants centralized sparse input in triplet form for the simplest C setup.
    std::vector<int> irn;
    std::vector<int> jcn;
    std::vector<double> a;
#endif
};

BEGIN_FECORE_CLASS(MumpsSparseSolver, LinearSolver)
    ADD_PARAMETER(m->printLevel, "print_level");
    ADD_PARAMETER(m->matrixMode, "matrix_mode");
END_FECORE_CLASS()

MumpsSparseSolver::MumpsSparseSolver(FEModel* fem) : LinearSolver(fem)
{
    m = new Impl();
    m->fem = fem;
}

MumpsSparseSolver::~MumpsSparseSolver()
{
    Destroy();
    delete m;
}

void MumpsSparseSolver::SetPrintLevel(int n)
{
    m->printLevel = n;
}

void MumpsSparseSolver::SetMatrixMode(int mode)
{
    m->matrixMode = mode;
}

SparseMatrix* MumpsSparseSolver::CreateSparseMatrix(Matrix_Type ntype)
{
    switch (ntype)
    {
        case REAL_SYMMETRIC: { m->isSymmetric = true; m->A = new CompactSymmMatrix(1); break; }
        case REAL_UNSYMMETRIC: { m->isSymmetric = false; m->A = new CRSSparseMatrix(1); break; }
        case REAL_SYMM_STRUCTURE: { m->isSymmetric = false; m->A = new CRSSparseMatrix(1); break; }
    default:
        assert(false);
        return nullptr;
    }
}

bool MumpsSparseSolver::SetSparseMatrix(SparseMatrix* A)
{
    if (m->isFactored) Destroy();

    m->A = dynamic_cast<CompactMatrix*>(A);
    if (m->A == nullptr) return false;

    m->isSymmetric = (m->A->isSymmetric() != 0);
    return true;
}

#ifdef HAVE_MUMPS
static int DetermineSym(MumpsSparseSolver::Impl* m)
{
    if (m->matrixMode == MumpsSparseSolver::MM_SYMMETRIC_POSDEF) return 1;
    if (m->matrixMode == MumpsSparseSolver::MM_SYMMETRIC_INDEF ) return 2;
    if (m->matrixMode == MumpsSparseSolver::MM_UNSYMMETRIC     ) return 0;

    return (m->isSymmetric ? 2 : 0);
}

static void ConvertCSCToTriplet(CompactMatrix* A,
                                std::vector<int>& irn,
                                std::vector<int>& jcn,
                                std::vector<double>& aval)
{
    const int n = A->Rows();
    const int nnz = A->NonZeroes();
    const int base = A->Offset();

    const int* colp = A->Pointers();
    const int* rowi = A->Indices();
    const double* ax = A->Values();

    irn.resize(nnz);
    jcn.resize(nnz);
    aval.resize(nnz);

    int k = 0;
    for (int col = 0; col < n; ++col)
    {
        int p0 = colp[col] - base;
        int p1 = colp[col + 1] - base;
        for (int p = p0; p < p1; ++p, ++k)
        {
            irn[k] = (rowi[p] - base) + 1; // MUMPS C interface uses 1-based sparse indices
            jcn[k] = col + 1;
            aval[k] = ax[p];
        }
    }
}
#endif

bool MumpsSparseSolver::PreProcess()
{
    assert(m->A);
    m->n = m->A->Rows();
    m->nnz = m->A->NonZeroes();

#ifndef HAVE_MUMPS
    feLogError("MUMPS support is not enabled in this build.");
    return false;
#else
    if (m->initialized)
    {
        m->id.job = JOB_END;
        dmumps_c(&m->id);
        m->initialized = false;
    }

    std::memset(&m->id, 0, sizeof(m->id));
    m->id.par = 1;     // host participates
    m->id.sym = DetermineSym(m);
    m->id.comm_fortran = -987654; // standard sequential/libmpiseq convention in many MUMPS C examples

    m->id.job = JOB_INIT;
    dmumps_c(&m->id);
    m->initialized = true;

    // Quiet by default unless debugging
    m->id.icntl[0] = -1;
    m->id.icntl[1] = -1;
    m->id.icntl[2] = -1;
    m->id.icntl[3] = (m->printLevel > 0 ? 6 : 0);

    ConvertCSCToTriplet(m->A, m->irn, m->jcn, m->a);

    m->id.n   = m->n;
    m->id.nz  = m->nnz;
    m->id.irn = m->irn.data();
    m->id.jcn = m->jcn.data();
    m->id.a   = m->a.data();

    m->id.job = 1; // analysis
    dmumps_c(&m->id);

    if (m->id.info[0] < 0)
    {
        feLogError("MUMPS analysis failed: INFOG(1)=%d INFOG(2)=%d",
            m->id.infog[0], m->id.infog[1]);
        return false;
    }

    return LinearSolver::PreProcess();
#endif
}

bool MumpsSparseSolver::Factor()
{
#ifndef HAVE_MUMPS
    feLogError("MUMPS support is not enabled in this build.");
    return false;
#else
    m->id.job = 2; // factorization
    dmumps_c(&m->id);

    if (m->id.info[0] < 0)
    {
        feLogError("MUMPS factorization failed: INFOG(1)=%d INFOG(2)=%d",
            m->id.infog[0], m->id.infog[1]);
        return false;
    }

    m->isFactored = true;
    return true;
#endif
}

bool MumpsSparseSolver::BackSolve(double* x, double* b)
{
#ifndef HAVE_MUMPS
    feLogError("MUMPS support is not enabled in this build.");
    return false;
#else
    if (!m->initialized || !m->isFactored) return false;

    if (x != b) std::memcpy(x, b, sizeof(double) * m->n);

    m->id.nrhs = 1;
    m->id.lrhs = m->n;
    m->id.rhs  = x;

    m->id.job = 3; // solve with existing factors
    dmumps_c(&m->id);

    if (m->id.info[0] < 0)
    {
        feLogError("MUMPS solve failed: INFOG(1)=%d INFOG(2)=%d",
            m->id.infog[0], m->id.infog[1]);
        return false;
    }

    return true;
#endif
}

void MumpsSparseSolver::Destroy()
{
#ifdef HAVE_MUMPS
    if (m->initialized)
    {
        m->id.job = JOB_END;
        dmumps_c(&m->id);
        m->initialized = false;
    }

    m->irn.clear();
    m->jcn.clear();
    m->a.clear();
#endif
    m->isFactored = false;
}

void MumpsSparseSolver::SetICNTL(int idx, int val)
{
#ifdef HAVE_MUMPS
    if (m->initialized && idx >= 1 && idx <= 60)
        m->id.icntl[idx - 1] = val;
#else
    (void)idx; (void)val;
#endif
}

void MumpsSparseSolver::SetCNTL(int idx, double val)
{
#ifdef HAVE_MUMPS
    if (m->initialized && idx >= 1 && idx <= 15)
        m->id.cntl[idx - 1] = val;
#else
    (void)idx; (void)val;
#endif
}
