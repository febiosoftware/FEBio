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
#include "LapackSolver.h"
#include "MatrixTools.h"
#include <FECore/log.h>

#ifdef __APPLE__

//////////////////////////////////////////////////////////////
// LapackSolver
//////////////////////////////////////////////////////////////

BEGIN_FECORE_CLASS(LapackSolver, LinearSolver)
ADD_PARAMETER(m_print_cn, "print_condition_number");
ADD_PARAMETER(m_iparm3  , "precondition");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
LapackSolver::LapackSolver(FEModel* fem) : LinearSolver(fem), m_pA(0)
{
    m_print_cn = false;
    m_mtype = -2;
    m_iparm3 = false;
    m_isFactored = false;
    
}

//-----------------------------------------------------------------------------
LapackSolver::~LapackSolver()
{
    Destroy();
}

//-----------------------------------------------------------------------------
void LapackSolver::PrintConditionNumber(bool b)
{
    m_print_cn = b;
}

//-----------------------------------------------------------------------------
void LapackSolver::UseIterativeFactorization(bool b)
{
    m_iparm3 = b;
}

//-----------------------------------------------------------------------------
SparseMatrix* LapackSolver::CreateSparseMatrix(Matrix_Type ntype)
{
    // allocate the correct matrix format depending on matrix symmetry type
    switch (ntype)
    {
        case REAL_SYMMETRIC     : m_mtype = -2; m_pA = new CompactSymmMatrix(0); break;
        case REAL_UNSYMMETRIC   : m_mtype = 11; m_pA = new CCSSparseMatrix(0); break;
        case REAL_SYMM_STRUCTURE: m_mtype =  1; m_pA = new CCSSparseMatrix(0); break;
        default:
            assert(false);
            m_pA = nullptr;
    }
    
    return m_pA;
}

//-----------------------------------------------------------------------------
bool LapackSolver::SetSparseMatrix(SparseMatrix* pA)
{
    if (m_pA && m_isFactored) Destroy();
    m_pA = dynamic_cast<CompactMatrix*>(pA);
    m_mtype = -2;
    if (dynamic_cast<CCSSparseMatrix*>(pA)) m_mtype = 11;
    return (m_pA != nullptr);
}

//-----------------------------------------------------------------------------
bool LapackSolver::PreProcess()
{
    assert(m_isFactored == false);
    
    m_n = m_pA->Rows();
    m_nnz = m_pA->NonZeroes();
    m_nrhs = 1;
    
    pointers.resize(m_nnz);
    for (int i=0; i<m_nnz; ++i) pointers[i] = m_pA->Pointers()[i];
    colS = &pointers[0];
    
    if (__builtin_available(macOS 10.13, *)) {
        SMS.columnCount = m_pA->Columns();
        SMS.rowCount = m_pA->Rows();
        SMS.rowIndices = m_pA->Indices();
        SMS.columnStarts = colS;
        SMS.blockSize = 1;
        if (m_mtype == 11) {
            // Use QR factorization for non-symmetric matrix
            SMS.attributes.kind = SparseOrdinary;
            ASS = SparseFactor(SparseFactorizationQR, SMS);
        }
        else {
            // Use LDLT factorization for symmetric matrix
            // (Cholesky factorization does not work)
            SMS.attributes.kind = SparseSymmetric;
            ASS = SparseFactor(SparseFactorizationLDLTTPP, SMS);
        }
    }
    else {
        // Fallback on earlier versions
        fprintf(stderr, "\nERROR during preprocessing: ");
        fprintf(stderr, "\nLAPACK solver not available for macOS earlier than 10.13!");
        return false;
    }
    
    return LinearSolver::PreProcess();
}

//-----------------------------------------------------------------------------
bool LapackSolver::Factor()
{
    // make sure we have work to do
    if (m_pA->Rows() == 0) return true;
    
    // ------------------------------------------------------------------------------
    // This step does the factorization
    // ------------------------------------------------------------------------------

    if (__builtin_available(macOS 10.13, *)) {
        // create the sparse matrix
        A.structure = SMS;
        A.data = m_pA->Values();
        
        ASF = SparseFactor(ASS,A);
        if (ASF.status != SparseStatusOK) {
            fprintf(stderr, "\nERROR during factorization:");
            switch (ASF.status) {
                case SparseFactorizationFailed:
                    fprintf(stderr, "\nFactorization failed!");
                    break;
                case SparseMatrixIsSingular:
                    fprintf(stderr, "\nSparse matrix is singular!");
                    break;
                case SparseInternalError:
                    fprintf(stderr, "\nSolver called internal error!");
                    break;
                case SparseParameterError:
                    fprintf(stderr, "\nSolver called parameter error!");
                    break;
                default:
                    fprintf(stderr, "\nUnknown error!");
            }
            return false;
        }
    }
    else {
        // Fallback on earlier versions
        fprintf(stderr, "\nERROR during factorization: ");
        fprintf(stderr, "\nLAPACK solver not available for macOS earlier than 10.13!");
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
bool LapackSolver::BackSolve(double* x, double* b)
{
    // make sure we have work to do
    if (m_pA->Rows() == 0) return true;
    if (m_isFactored == false) return true;
    
    DenseVector_Double X, B;
    X.count = B.count = m_n;
    X.data = x; B.data = b;
    
    // solve the system
    if (__builtin_available(macOS 10.13, *)) {
        SparseSolve(ASF,B,X);
    } else {
        // Fallback on earlier versions
        fprintf(stderr, "\nERROR during back solve: ");
        fprintf(stderr, "\nLAPACK solver not available for macOS earlier than 10.13!");
        return false;
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
double LapackSolver::condition_number()
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
        
        double normx = NumCore::infNorm(x);
        if (normx > normAi) normAi = normx;
        
        int pct = (100 * i) / (iters - 1);
        fprintf(stderr, "calculating condition number: %d%%\r", pct);
    }
    
    double c = normA*normAi;
    return c;
}

//-----------------------------------------------------------------------------
void LapackSolver::Destroy()
{
    if (m_pA && m_pA->Pointers() && m_isFactored)
    {
        SparseCleanup(ASS);
        SparseCleanup(ASF);
    }
    m_isFactored = false;
}
#else
BEGIN_FECORE_CLASS(LapackSolver, LinearSolver)
ADD_PARAMETER(m_print_cn, "print_condition_number");
ADD_PARAMETER(m_iparm3, "precondition");
END_FECORE_CLASS();

LapackSolver::LapackSolver(FEModel* fem) : LinearSolver(fem) {}
LapackSolver::~LapackSolver() {}
bool LapackSolver::PreProcess() { return false; }
bool LapackSolver::Factor() { return false; }
bool LapackSolver::BackSolve(double* x, double* y) { return false; }
void LapackSolver::Destroy() {}
SparseMatrix* LapackSolver::CreateSparseMatrix(Matrix_Type ntype) { return nullptr; }
bool LapackSolver::SetSparseMatrix(SparseMatrix* pA) { return false; }
void LapackSolver::PrintConditionNumber(bool b) {}
double LapackSolver::condition_number() { return 0; }
void LapackSolver::UseIterativeFactorization(bool b) {}
#endif
