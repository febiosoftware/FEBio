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



#pragma once
#include <FECore/LinearSolver.h>

class BoomerAMGSolver : public LinearSolver
{
	class Implementation;

public:
	enum InterpolationType {
		INTERPOL_OP_0 = 0,			// classicial modified interpolation
		INTERPOL_OP_1 = 1,			// LS interpolation (for use with GSMG)
		INTERPOL_OP_2 = 2,			// classical modified interpolation for hyperbolic PDEs
		INTERPOL_OP_3 = 3,			// direct interpolation (with separation of weights)
		INTERPOL_OP_4 = 4,			// multipass interpolation
		INTERPOL_OP_5 = 5,			// multipass interpolation (with separation of weights)
		INTERPOL_OP_6 = 6,			// extended+i interpolation	
		INTERPOL_OP_7 = 7,			// extended+i (if no common C neighbor) interpolation
		INTERPOL_OP_8 = 8,			// standard interpolation
		INTERPOL_OP_9 = 9,			// standard interpolation (with separation of weights)
		INTERPOL_OP_10 = 10,		// classical block interpolation (for use with nodal systems version only)
		INTERPOL_OP_11 = 11,		// like 10, with diagonalized blocks
		INTERPOL_OP_12 = 12,		// FF interpolation
		INTERPOL_OP_13 = 13,		// FF1 interpolation
		INTERPOL_OP_14 = 14,		// extended interpolation
	};

	enum RelaxationType {
		JACOBI = 0,							// Jacobi
		GAUSS_SEIDEL_SEQ = 1,				// Gauss-Seidel, sequential (very slow!)
		GAUSS_SEIDEL_PAR = 2,				// Gauss-Seidel, interior points in parallel, boundary sequential (slow!)
		SOR_FORWARD = 3,					// hybrid Gauss-Seidel or SOR, forward solve
		SOR_BACKWARD = 4,					// hybrid Gauss-Seidel or SOR, backward solve
		CHAOTIC_GUASS_SEIDEL = 5,			// hybrid chaotic Gauss-Seidel (works only with OpenMP)
		SYMMETRIC_GAUSS_SEIDEL = 6,			// hybrid symmetric Gauss-Seidel or SSOR
		SYMMETRIC_GAUSS_SEIDEL_L1 = 8,		// l1-scaled hybrid symmetric Gauss-Seidel
		GAUSSIAN_ELIMINATION = 9,			// Gaussian elimination (only on coarsest level)
		CG = 15,							// CG (warning - not a fixed smoother - may require RGMRES)
		CHEBYSHEV  = 16,					// Chebyshev
		FCF_JACOBI = 17,					// FCF-Jacobi
		JACOBI_L1  = 18						// l1-scaled Jacobi
	};

	enum CoarsenType {
		COARSE_OP_0 = 0,			// CLJP-coarsening (a parallel coarsening algorithm using independent sets)
		COARSE_OP_1 = 1,			// classical Ruge-Stueben coarsening, no boundary treatment (not recommended!)
		COARSE_OP_3 = 3,			// classical Ruge-Stueben coarsening, with boundary treatment
		COARSE_OP_6 = 6,			// Falgout coarsening
		COARSE_OP_7 = 7,			// CLJP-coarsening (using fix random vector, for debugging purposes only)
		COARSE_OP_8 = 8,			// PMIS-coarsening (might lead to slower convergence)
		COARSE_OP_9 = 9,			// PMIS-coarsening (using a fixed random vector ,for debugging purposes only)
		COARSE_OP_10 = 10,			// HMIS-coarsening
		COARSE_OP_11 = 11,			// one-pass Ruge-Stueben coarsening (not recommended!)
		COARSE_OP_21 = 21,			// CGC coarsening
		COARSE_OP_22 = 22,			// CGC-E coarsening
	};
    
    enum AggInterpType {
        AGG_IT_0 = 0,
        AGG_IT_5 = 5,
        AGG_IT_7 = 7,
    };

public:
	BoomerAMGSolver(FEModel* fem);
	~BoomerAMGSolver();

public:
	void SetPrintLevel(int printLevel) override;
	void SetMaxIterations(int maxIter);
	void SetConvergenceTolerance(double tol);
	void SetMaxLevels(int levels);
	void SetCoarsenType(int coarsenType);
	void SetUseNumFunctions(bool b);
    void SetRelaxType(int rlxtyp);
    void SetInterpType(int inptyp);
    void SetStrongThreshold(double thresh);
    void SetPMaxElmts(int pmax);
    void SetNumSweeps(int nswp);
    void SetAggInterpType(int aggit);
    void SetAggNumLevels(int anlv);
	void SetNodal(int nodal);
	void SetJacobiPC(bool b);
	void SetFailOnMaxIterations(bool b);
	bool GetJacobiPC();

public:
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
	bool SetSparseMatrix(SparseMatrix* pA) override;

	bool PreProcess() override;
	bool Factor() override;
	bool BackSolve(double* x, double* y) override;
	void Destroy() override;

private:
	Implementation*	imp;

	DECLARE_FECORE_CLASS();
};
