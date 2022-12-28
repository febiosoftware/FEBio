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
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif

// This sparse linear solver uses the Accelerate framework on MAC. 
class AccelerateSparseSolver : public LinearSolver
{
    class   Implementation;
    
public:
    enum FactorizationType {
        FTSparseFactorizationCholesky = 0,
        FTSparseFactorizationLDLT = 1,
        FTSparseFactorizationLDLTUnpivoted = 2,
        FTSparseFactorizationLDLTSBK = 3,
        FTSparseFactorizationLDLTTPP = 4,
        FTSparseFactorizationQR = 5,
        FTSparseFactorizationCholeskyAtA = 6,
    };
    
    enum OrderMethod {
        OMSparseOrderAMD = 0,
        OMSparseOrderMetis = 1,
        OMSparseOrderCOLAMD = 2,
    };

    enum IterativeMethod {
        ITSparseConjugateGradient = 0,
        ITSparseGMRES = 1,
        ITSparseDQGMRES = 2,
        ITSparseFGMRES = 3,
        ITSparseLSMR = 4,
    };
    
public:
    AccelerateSparseSolver(FEModel* fem);
    ~AccelerateSparseSolver();
    
public:
    void PrintConditionNumber(bool b);
    void UseIterativeFactorization(bool b);
    void SetFactorizationType(int ftype);
    void SetOrderMethod(int order);
    void SetPrintLevel(int printlevel) override;

    double condition_number();
    static void MyReportError(const char* message);
    static void MyReportStatus(const char* message);

public:
    SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
    bool SetSparseMatrix(SparseMatrix* pA) override;
    
    bool PreProcess() override;
    bool Factor() override;
    bool BackSolve(double* x, double* y) override;
    void Destroy() override;
    bool IsIterative() const override;
    
private:
    Implementation*    imp;
        
    DECLARE_FECORE_CLASS();
};
