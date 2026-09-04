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
#include <FECore/CompactMatrix.h>
#include "numcore_api.h"

class NUMCORE_API MumpsSparseSolver : public LinearSolver
{
public:
    enum MatrixMode {
        MM_AUTO = 0,
        MM_SYMMETRIC_POSDEF = 1,
        MM_SYMMETRIC_INDEF  = 2,
        MM_UNSYMMETRIC      = 3
    };

public:
    MumpsSparseSolver(FEModel* fem);
    ~MumpsSparseSolver() override;

    void SetPrintLevel(int n) override;
    void SetMatrixMode(int mode);
    void SetICNTL(int idx, int val);
    void SetCNTL(int idx, double val);

public:
    SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
    bool SetSparseMatrix(SparseMatrix* A) override;
    bool PreProcess() override;
    bool Factor() override;
    bool BackSolve(double* x, double* b) override;
    void Destroy() override;
    bool IsIterative() const override { return false; }
    
public:
    class Impl;
    Impl* m;
    
    DECLARE_FECORE_CLASS();
};
