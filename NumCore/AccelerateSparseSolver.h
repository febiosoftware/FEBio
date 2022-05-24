/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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



#pragma once
#include <FECore/LinearSolver.h>
#include "CompactUnSymmMatrix.h"
#include "CompactSymmMatrix.h"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif

// This sparse linear solver uses the Accelerate framework on MAC. 
class AccelerateSparseSolver : public LinearSolver
{
public:
    AccelerateSparseSolver(FEModel* fem);
    ~AccelerateSparseSolver();
    bool PreProcess() override;
    bool Factor() override;
    bool BackSolve(double* x, double* y) override;
    void Destroy() override;
    
    SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
    bool SetSparseMatrix(SparseMatrix* pA) override;
    
    void PrintConditionNumber(bool b);
    
    double condition_number();
    
    void UseIterativeFactorization(bool b);
    
protected:
    
    CompactMatrix*    m_pA;
    int                m_mtype; // matrix type
    
#ifdef __APPLE__
    SparseMatrixStructure SMS;
    SparseMatrix_Double A;
    SparseOpaqueSymbolicFactorization ASS;
    SparseOpaqueFactorization_Double ASF;
    long* colS;
    std::vector<long> pointers;
#endif
    
    // solver control parameters
    
    bool m_iparm3;    // use direct-iterative method
    
    // Matrix data
    int m_n, m_nnz, m_nrhs;
    
    bool    m_print_cn;    // estimate and print the condition number
    
    bool    m_isFactored;
    
    DECLARE_FECORE_CLASS();
};
