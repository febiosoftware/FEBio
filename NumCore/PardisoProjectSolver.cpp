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
#include <stdio.h>
#include <stdlib.h>
#include "PardisoProjectSolver.h"
#include "MatrixTools.h"
#include <FECore/log.h>

//! This implementation of the Pardiso solver is for the version
//! available in the Intel MKL.


#ifdef PARDISODL

#include <thread>

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                             double *, int    *,    int *, int *,   int *, int *,
                             int *, double *, double *, int *, double *);

//-----------------------------------------------------------------------------
// print pardiso error message
void print_err_pdl(int nerror)
{
    switch (-nerror)
    {
        case 1: fprintf(stderr, "Inconsistent input\n"); break;
        case 2: fprintf(stderr, "Not enough memory\n"); break;
        case 3: fprintf(stderr, "Reordering problem\n"); break;
        case 4: fprintf(stderr, "Zero pivot, numerical fact. or iterative refinement problem\n"); break;
        case 5: fprintf(stderr, "Unclassified (internal) error\n"); break;
        case 6: fprintf(stderr, "Preordering failed\n"); break;
        case 7: fprintf(stderr, "Diagonal matrix problem\n"); break;
        case 8: fprintf(stderr, "32-bit integer overflow problem\n"); break;
        default:
            fprintf(stderr, " Unknown\n");
    }
}

//////////////////////////////////////////////////////////////
// PardisoProjectSolver
//////////////////////////////////////////////////////////////

BEGIN_FECORE_CLASS(PardisoProjectSolver, LinearSolver)
ADD_PARAMETER(m_print_cn, "print_condition_number");
ADD_PARAMETER(m_iparm3  , "precondition");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
PardisoProjectSolver::PardisoProjectSolver(FEModel* fem) : LinearSolver(fem), m_pA(0)
{
    m_print_cn = false;
    m_mtype = -2;
    m_iparm3 = false;
    m_isFactored = false;
    
    /* If both PARDISO AND PARDISODL are defined, print a warning */
#ifdef PARDISO
    fprintf(stderr, "WARNING: The pardiso-project version of the Pardiso solver is being used\n\n");
    exit(1);
#endif
}

//-----------------------------------------------------------------------------
PardisoProjectSolver::~PardisoProjectSolver()
{
}

//-----------------------------------------------------------------------------
void PardisoProjectSolver::PrintConditionNumber(bool b)
{
    m_print_cn = b;
}

//-----------------------------------------------------------------------------
void PardisoProjectSolver::UseIterativeFactorization(bool b)
{
    m_iparm3 = b;
}

//-----------------------------------------------------------------------------
SparseMatrix* PardisoProjectSolver::CreateSparseMatrix(Matrix_Type ntype)
{
    // allocate the correct matrix format depending on matrix symmetry type
    switch (ntype)
    {
        case REAL_SYMMETRIC     : m_mtype = -2; m_pA = new CompactSymmMatrix(1); break;
        case REAL_UNSYMMETRIC   : m_mtype = 11; m_pA = new CRSSparseMatrix(1); break;
        case REAL_SYMM_STRUCTURE: m_mtype =  1; m_pA = new CRSSparseMatrix(1); break;
        default:
            assert(false);
            m_pA = nullptr;
    }
    
    return m_pA;
}

//-----------------------------------------------------------------------------
bool PardisoProjectSolver::SetSparseMatrix(SparseMatrix* pA)
{
    if (m_pA && m_isFactored) Destroy();
    m_pA = dynamic_cast<CompactMatrix*>(pA);
    m_mtype = -2;
    if (dynamic_cast<CRSSparseMatrix*>(pA)) m_mtype = 11;
    return (m_pA != nullptr);
}

//-----------------------------------------------------------------------------
bool PardisoProjectSolver::PreProcess()
{
    m_iparm[0] = 0; /* Use default values for parameters */
//    for (int i=0; i<64; ++i) m_pt[i] = nullptr;
    
    //fprintf(stderr, "In PreProcess\n");
    assert(m_isFactored == false);
    int error = 0;
    int solver = 0; /* use sparse direct solver */
    pardisoinit(m_pt, &m_mtype, &solver, m_iparm, m_dparm, &error);
    
    m_n = m_pA->Rows();
    m_nnz = m_pA->NonZeroes();
    m_nrhs = 1;
    
    /* Numbers of processors, value of OMP_NUM_THREADS (if available) or
     total number of processors on machine */
    const auto processor_count = std::thread::hardware_concurrency();
    int num_procs = processor_count;
    char* var = getenv("OMP_NUM_THREADS");
    if (var != NULL)
        sscanf( var, "%d", &num_procs );
    feLog("Number of processors: %d\n",num_procs);
    if (num_procs == 0) {
        feLogError("Specify number of processors in environmental variable OMP_NUM_THREADS!");
        return false;
    }
    m_iparm[2]  = num_procs;

    m_maxfct = 1;    /* Maximum number of numerical factorizations */
    m_mnum = 1;    /* Which factorization to use */
    
    m_msglvl = 0;    /* 0 Suppress printing, 1 Print statistical information */
    
    return LinearSolver::PreProcess();
}

//-----------------------------------------------------------------------------
bool PardisoProjectSolver::Factor()
{
    // make sure we have work to do
    if (m_pA->Rows() == 0) return true;
    
    // ------------------------------------------------------------------------------
    // Reordering and Symbolic Factorization.  This step also allocates all memory
    // that is necessary for the factorization.
    // ------------------------------------------------------------------------------
    
    int phase = 11;
    
    int error = 0;
    pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, m_pA->Values(), m_pA->Pointers(), m_pA->Indices(),
            NULL, &m_nrhs, m_iparm, &m_msglvl, NULL, NULL, &error, m_dparm);
    
    if (error)
    {
        fprintf(stderr, "\nERROR during symbolic factorization: ");
        print_err_pdl(error);
        exit(2);
    }
    
    // ------------------------------------------------------------------------------
    // This step does the factorization
    // ------------------------------------------------------------------------------
    
    phase = 22;
    
    m_iparm[3] = (m_iparm3 ? 61 : 0);
    error = 0;
    pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, m_pA->Values(), m_pA->Pointers(), m_pA->Indices(),
            NULL, &m_nrhs, m_iparm, &m_msglvl, NULL, NULL, &error, m_dparm);
    
    if (error)
    {
        fprintf(stderr, "\nERROR during factorization: ");
        print_err_pdl(error);
        return false;
    }
    
    // calculate and print the condition number
    if (m_print_cn)
    {
        double c = condition_number();
        feLog("\tcondition number (est.) ................... : %lg\n\n", c);
    }
    
    m_isFactored = true;
    
    return true;
}

//-----------------------------------------------------------------------------
bool PardisoProjectSolver::BackSolve(double* x, double* b)
{
    // make sure we have work to do
    if (m_pA->Rows() == 0) return true;
    
    int phase = 33;
    
    m_iparm[7] = 1;    /* Maximum number of iterative refinement steps */
    
    int error = 0;
    pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, m_pA->Values(), m_pA->Pointers(), m_pA->Indices(),
            NULL, &m_nrhs, m_iparm, &m_msglvl, b, x, &error, m_dparm);
    
    if (error)
    {
        fprintf(stderr, "\nERROR during solution: ");
        print_err_pdl(error);
        exit(3);
    }
    
    // update stats
    UpdateStats(1);
    
    return true;
}

//-----------------------------------------------------------------------------
// This algorithm (naively) estimates the condition number. It is based on the observation that
// for a linear system of equations A.x = b, the following holds
// || A^-1 || >= ||x||.||b||
// Thus the condition number can be estimated by
// c = ||A||.||A^-1|| >= ||A|| . ||x|| / ||b||
// This algorithm tries for some random b vectors with norm ||b||=1 to maxize the ||x||.
// The returned value will be an underestimate of the condition number
double PardisoProjectSolver::condition_number()
{
    // This assumes that the factorization is already done!
    int N = m_pA->Rows();
    
    // get the norm of the matrix
    double normA = m_pA->infNorm();
    
    // estimate the norm of the inverse of A
    double normAi = 0.0;
    
    // choose max iterations
    int iters = (N < 50 ? N : 50);
    
    vector<double> b(N, 0), x(N, 0);
    for (int i = 0; i < iters; ++i)
    {
        // create a random vector
        NumCore::randomVector(b, -1.0, 1.0);
        for (int j = 0; j < N; ++j) b[j] = (b[j] >= 0.0 ? 1.0 : -1.0);
        
        // calculate solution
        BackSolve(&x[0], &b[0]);
        
        double normb = NumCore::infNorm(b);
        double normx = NumCore::infNorm(x);
        if (normx > normAi) normAi = normx;
        
        int pct = (100 * i) / (iters - 1);
        fprintf(stderr, "calculating condition number: %d%%\r", pct);
    }
    
    double c = normA*normAi;
    return c;
}

//-----------------------------------------------------------------------------
void PardisoProjectSolver::Destroy()
{
    int phase = -1;
    
    int error = 0;
    
    if (m_pA && m_pA->Pointers() && m_isFactored)
    {
        pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, NULL, m_pA->Pointers(), m_pA->Indices(),
                NULL, &m_nrhs, m_iparm, &m_msglvl, NULL, NULL, &error, m_dparm);
    }
    m_isFactored = false;
}

#else
BEGIN_FECORE_CLASS(PardisoProjectSolver, LinearSolver)
	ADD_PARAMETER(m_print_cn, "print_condition_number");
	ADD_PARAMETER(m_iparm3, "precondition");
END_FECORE_CLASS();

PardisoProjectSolver::PardisoProjectSolver(FEModel* fem) : LinearSolver(fem) {}
PardisoProjectSolver::~PardisoProjectSolver() {}
bool PardisoProjectSolver::PreProcess() { return false; }
bool PardisoProjectSolver::Factor() { return false; }
bool PardisoProjectSolver::BackSolve(double* x, double* y) { return false; }
void PardisoProjectSolver::Destroy() {}
SparseMatrix* PardisoProjectSolver::CreateSparseMatrix(Matrix_Type ntype) { return nullptr; }
bool PardisoProjectSolver::SetSparseMatrix(SparseMatrix* pA) { return false; }
void PardisoProjectSolver::PrintConditionNumber(bool b) {}
double PardisoProjectSolver::condition_number() { return 0; }
void PardisoProjectSolver::UseIterativeFactorization(bool b) {}
#endif
