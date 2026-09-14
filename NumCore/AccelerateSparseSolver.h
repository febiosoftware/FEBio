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
    //! Factorization types.
    //
    //  IMPORTANT: these values are the FEBio-facing encoding. They are NOT
    //  Accelerate's SparseFactorization_t values (where QR = 40 and
    //  CholeskyAtA = 41); toSparseFactorization() in the .cpp translates.
    //
    //  They match the numbering that FEBio input files have always used, so
    //  existing files keep working. Input files should prefer the names in the
    //  comments below - those are resolved by ResolveParameters(), which reads
    //  the option as a string and accepts either the name or the integer.
    enum FactorizationType {
        FTSparseFactorizationCholesky       = 0,   // "cholesky"
        FTSparseFactorizationLDLT           = 1,   // "ldlt"
        FTSparseFactorizationLDLTUnpivoted  = 2,   // "ldlt_unpivoted"
        FTSparseFactorizationLDLTSBK        = 3,   // "ldlt_sbk"
        FTSparseFactorizationLDLTTPP        = 4,   // "ldlt_tpp"
        FTSparseFactorizationQR             = 5,   // "qr"
        FTSparseFactorizationCholeskyAtA    = 6,   // "cholesky_ata"
        FTSparseFactorizationDefault        = 7,   // "default"
        // LU requires macOS 15.5 / iOS 18.5 or later, both to build and to run.
        FTSparseFactorizationLU             = 8,   // "lu"            (default LU, currently TPP)
        FTSparseFactorizationLUUnpivoted    = 9,   // "lu_unpivoted"
        FTSparseFactorizationLUSPP          = 10,  // "lu_spp"        (pivoting within supernodes)
        FTSparseFactorizationLUTPP          = 11,  // "lu_tpp"        (threshold partial pivoting)
    };

    //! Scaling applied by Accelerate itself during the numeric factorization.
    //
    //  This is independent of, and composes with, the partition scaling below.
    //  Accelerate's version is a global equilibration; it knows nothing about
    //  the velocity/dilatation block structure. For LU the framework default is
    //  no scaling, so "equilibriation_inf" is worth trying on ill-conditioned
    //  systems before reaching for anything custom.
    enum FactorScaling {
        FSDefault           = 0,   // "default"             (LDLT: equilibriation; Cholesky/LU: none)
        FSNone              = 1,   // "none"
        FSEquilibriationInf = 2,   // "equilibriation_inf"
    };

    enum OrderMethod {
        OMSparseOrderAMD     = 0,   // "amd"
        OMSparseOrderMetis   = 1,   // "metis"
        OMSparseOrderCOLAMD  = 2,   // "colamd"
        OMSparseOrderDefault = 3,   // "default"
    };

    enum IterativeMethod {
        ITSparseConjugateGradient = 0,   // "cg"
        ITSparseGMRES             = 1,   // "gmres"
        ITSparseDQGMRES           = 2,   // "dqgmres"
        ITSparseFGMRES            = 3,   // "fgmres"
        ITSparseLSMR              = 4,   // "lsmr"
    };

    //! Preconditioner used by the iterative methods.
    enum PreconditionerType {
        PCNone         = 0,   // "none"
        PCDiagonal     = 1,   // "diagonal"
        PCDiagScaling  = 2,   // "diag_scaling"
    };

    //! How the linear system is solved.
    //
    //  SSFactorPreconditioned keeps a direct factorization around and uses it
    //  as the preconditioner for GMRES rather than as the solver. The stored
    //  factorization is allowed to go stale: an out-of-date factorization is
    //  usually a very good preconditioner even when it is a poor solver, so
    //  GMRES converges in a handful of iterations and the matrix only has to be
    //  refactored occasionally. This is the linear-algebra analogue of what
    //  Broyden's method already does at the nonlinear level.
    enum SolveStrategy {
        SSDirect               = 0,   // "direct"
        SSIterative            = 1,   // "iterative"
        SSFactorPreconditioned = 2,   // "factor_preconditioned"
    };

    //! Partition-aware diagonal scaling of the linear system.
    //
    //  Kept as an enum rather than a bool so that call sites read as
    //  setPartitionScaling(PartitionScalingOn) instead of (true), and so that
    //  additional strategies can be added without an API break.
    enum PartitionScaling {
        PartitionScalingOff = 0,   // "off"
        PartitionScalingOn  = 1,   // "on"
    };

public:
    AccelerateSparseSolver(FEModel* fem);
    ~AccelerateSparseSolver();

public:
    void PrintConditionNumber(bool b);
    void UseIterativeFactorization(bool b);
    void SetFactorizationType(int ftype);
    void SetOrderMethod(int order);
    void SetIterativeMethod(int method);
    void SetPartitionScaling(PartitionScaling s);
    void SetSolveStrategy(int strategy);
    void SetPrintLevel(int printlevel) override;

    double condition_number();
    static void MyReportError(const char* message);
    static void MyReportStatus(const char* message);

#ifdef __APPLE__
    //! Preconditioner callback for SSFactorPreconditioned. Computes Y = P X
    //! where P is the stored (possibly stale) factorization of A.
    static void ApplyStoredFactor(void* mem, enum CBLAS_TRANSPOSE trans,
                                  DenseMatrix_Double X, DenseMatrix_Double Y);
#endif

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

    //! Build the partition scaling vector D and the scaled value array
    //! A' = D A D. Called once per Factor().
    bool ComputeAndApplyScaling();

    //! Perform the numeric factorization of the current matrix. Used by
    //! Factor() and by the stale-preconditioner retry in BackSolve().
    bool FactorNow();

    //! Resolve the string-valued options ("qr", "gmres", "on", ...) into the
    //! integer members. Called at the top of PreProcess().
    bool ResolveParameters();

    DECLARE_FECORE_CLASS();
};
