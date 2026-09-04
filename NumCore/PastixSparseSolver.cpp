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
#include "PastixSparseSolver.h"
#include "MatrixTools.h"
#include <FECore/log.h>
#include <cstring>
#include <iostream>

PastixSparseSolver::PastixSparseSolver(FEModel* fem)
    : LinearSolver(fem)
{
}

PastixSparseSolver::~PastixSparseSolver()
{
    Destroy();
}

SparseMatrix* PastixSparseSolver::CreateSparseMatrix(Matrix_Type ntype)
{
    switch (ntype)
    {
        case REAL_SYMMETRIC:    m_mtype = SpmSymmetric; m_pA = new CompactSymmMatrix(1); break;
        case REAL_UNSYMMETRIC:
        case REAL_SYMM_STRUCTURE: m_mtype = SpmGeneral; m_pA = new CRSSparseMatrix(1);
            break;
        default:
            m_pA = new CRSSparseMatrix(1);
            break;
    }

    return m_pA;
}

bool PastixSparseSolver::SetSparseMatrix(SparseMatrix* pA)
{
    if (m_pA && m_factored) Destroy();
    m_pA = dynamic_cast<CompactMatrix*>(pA);
    m_mtype = SpmSymmetric;
    if (dynamic_cast<CRSSparseMatrix*>(pA)) m_mtype = SpmGeneral;
    return (m_pA != nullptr);
}

bool PastixSparseSolver::PreProcess()
{
    if (m_pA == nullptr) return false;

    m_n = m_pA->Rows();
    m_nrhs = 1;

    pastixInitParam(m_iparm, m_dparm);

    m_iparm[IPARM_MODIFY_PARAMETER] = 1;
    m_iparm[IPARM_THREAD_NBR]       = 1;
    m_iparm[IPARM_VERBOSE]          = PastixVerboseNot;
    m_iparm[IPARM_FLOAT]            = PastixDouble;
    m_iparm[IPARM_MTX_TYPE]         = m_mtype;
    m_iparm[IPARM_FACTORIZATION]    = (m_mtype == SpmSymmetric) ? PastixFactLDLT : PastixFactLU;
    m_iparm[IPARM_REFINEMENT]       = PastixRefineGMRES;
    m_iparm[IPARM_ITERMAX]          = 100;
    m_iparm[IPARM_MC64]             = 1;
    m_iparm[IPARM_ORDERING]         = PastixOrderScotch;
    m_iparm[35]          = 2;    // enable internal two-sided scaling/equilibration

    pastixInit(&m_pastix, 0, m_iparm, m_dparm);

    m_initialized = true;
    m_factored = false;

    return LinearSolver::PreProcess();
}

bool PastixSparseSolver::CopyPatternToPastix()
{
    if (m_pA == nullptr) return false;

    const int n = m_pA->Rows();
    const int nnz = static_cast<int>(m_pA->NonZeroes());

    int* p = m_pA->Pointers();
    int* i = m_pA->Indices();

    if (p == nullptr || i == nullptr) return false;

    m_colptr.resize(n + 1);
    m_rowind.resize(nnz);
    m_perm.assign(n, 0);
    m_invp.assign(n, 0);

    const int offset = m_pA->Offset();

    for (int k = 0; k <= n; ++k)
        m_colptr[k] = static_cast<pastix_int_t>(p[k] - offset);

    for (int k = 0; k < nnz; ++k)
        m_rowind[k] = static_cast<pastix_int_t>(i[k] - offset);

    return true;
}

bool PastixSparseSolver::Factor()
{
    if (m_pA == nullptr) return false;
    if (m_pA->Rows() == 0) return true;

    if (!m_initialized) {
        if (!PreProcess()) return false;
    }

    if (!CopyPatternToPastix()) return false;

    std::vector<double> dummy_rhs(m_n, 0.0);

    m_iparm[IPARM_START_TASK] = PastixTaskOrdering;
    m_iparm[IPARM_END_TASK]   = PastixTaskNumfact;

    int ierr = pastix(
        &m_pastix,
        0,
        static_cast<pastix_int_t>(m_n),
        m_colptr.data(),
        m_rowind.data(),
        m_pA->Values(),
        m_perm.data(),
        m_invp.data(),
        dummy_rhs.data(),
        static_cast<pastix_int_t>(m_nrhs),
        m_iparm,
        m_dparm
    );

// For debugging
//    feLog("Number of refinement iterations: %d\n", m_iparm[IPARM_NBITER]);
    
    if (ierr != 0) {
        std::cerr << "PaStiX factorization error: " << ierr << std::endl;
        return false;
    }

    m_factored = true;
    return true;
}

bool PastixSparseSolver::BackSolve(double* x, double* b)
{
    if (m_pA == nullptr) return false;
    if (m_pA->Rows() == 0) return true;
    if (!m_factored) return false;

    // PaStiX overwrites RHS with the solution.
    std::vector<double> rhs(b, b + m_n);

    m_iparm[IPARM_START_TASK] = PastixTaskSolve;
    m_iparm[IPARM_END_TASK]   = PastixTaskRefine;

    int ierr = pastix(
        &m_pastix,
        0,
        static_cast<pastix_int_t>(m_n),
        m_colptr.data(),
        m_rowind.data(),
        m_pA->Values(),
        m_perm.data(),
        m_invp.data(),
        rhs.data(),
        static_cast<pastix_int_t>(m_nrhs),
        m_iparm,
        m_dparm
    );

    if (ierr != 0) {
        std::cerr << "PaStiX solve error: " << ierr << std::endl;
        return false;
    }

    std::memcpy(x, rhs.data(), sizeof(double) * m_n);

    UpdateStats(1);
    return true;
}

void PastixSparseSolver::Destroy()
{
    if (m_pastix != nullptr) {
        pastixFinalize(&m_pastix);
    }

    m_initialized = false;
    m_factored = false;

    m_colptr.clear();
    m_rowind.clear();
    m_perm.clear();
    m_invp.clear();
}
