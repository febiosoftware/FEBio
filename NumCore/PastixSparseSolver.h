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
#include <FECore/CompactUnSymmMatrix.h>
#include <FECore/CompactSymmMatrix.h>
#include "numcore_api.h"

#include <pastix.h>
#include "/usr/local/include/metis.h"
#include <vector>

class PastixSparseSolver : public LinearSolver
{
public:
    PastixSparseSolver(FEModel* fem);
    ~PastixSparseSolver();

    SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
    bool SetSparseMatrix(SparseMatrix* pA) override;

    bool PreProcess() override;
    bool Factor() override;
    bool BackSolve(double* x, double* b) override;
    void Destroy() override;

    double ConditionNumber() override { return 0.0; }

private:
    CompactMatrix* m_pA = nullptr;
    int m_mtype; // matrix type

    pastix_data_t* m_pastix = nullptr;
    pastix_int_t m_iparm[IPARM_SIZE];
    double m_dparm[DPARM_SIZE];

    std::vector<pastix_int_t> m_colptr;
    std::vector<pastix_int_t> m_rowind;
    std::vector<pastix_int_t> m_perm;
    std::vector<pastix_int_t> m_invp;

    int m_n = 0;
    int m_nrhs = 1;
    bool m_initialized = false;
    bool m_factored = false;

    std::vector<double> m_scale;        // diagonal scaling factors, size n
    std::vector<double> m_scaledValues; // scaled copy of matrix values, size nnz
    
    bool ComputeAndApplyScaling();
    
    bool CopyPatternToPastix();
};
