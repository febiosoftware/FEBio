/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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

class FGMRES_ILU0_Solver;

// Experimental solver class that can switch between direct (Pardiso) and 
// iterative (FGMRES\ILU0) solver.
// Currently, only supports non-symmetric matrices
class MixedLinearSolver : public LinearSolver
{
public:
	enum Strategy {
		DIRECT_SOLVER,
		ITERATIVE_SOLVER
	};

public:
	MixedLinearSolver(FEModel* fem);
	~MixedLinearSolver();
	bool PreProcess() override;
	bool Factor() override;
	bool BackSolve(double* x, double* y) override;
	void Destroy() override;

	void SetSolverStrategy(Strategy n);

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
	bool SetSparseMatrix(SparseMatrix* pA) override;

public: // properties for iterative solver
	void SetMaxIterations(int nmax);
	void SetPrintLevel(int n) override;
	void SetRelativeConvergence(double tol);
	void SetAbsoluteConvergence(double tol);

protected:
	LinearSolver* currentSolver() { return m_solver[m_strategy]; }

protected:
	int		m_strategy;	// 0 = direct, 1 = iterative
	CRSSparseMatrix*	m_A;
	LinearSolver*	m_solver[2];
	FGMRES_ILU0_Solver*	m_fgmres;
};
