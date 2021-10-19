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
#include "BlockMatrix.h"
#include "PardisoSolver.h"

//-----------------------------------------------------------------------------
// This class implements an iterative block solution strategy for solving linear 
// systems.

class BlockIterativeSolver : public IterativeLinearSolver
{
public:
	enum SOLUTION_METHOD {
		JACOBI,
		GAUSS_SEIDEL
	};

public:
	//! constructor
	BlockIterativeSolver(FEModel* fem);

	//! destructor
	~BlockIterativeSolver();

	//! Preprocess 
	bool PreProcess() override;

	//! Factor matrix
	bool Factor() override;

	//! Backsolve the linear system
	bool BackSolve(double* x, double* y) override;

	//! Clean up
	void Destroy() override;

	//! Create a sparse matrix
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	//! set the sparse matrix
	bool SetSparseMatrix(SparseMatrix* m) override;

	// return whether the iterative solver has a preconditioner or not
	bool HasPreconditioner() const override;

public:
	// Set the relative convergence tolerance
	void SetRelativeTolerance(double tol);

	// set the max nr of iterations
	void SetMaxIterations(int maxiter);

	// get the iteration count
	int GetIterations() const;

	// set fail on max iterations flag
	void SetFailMaxIters(bool b);

	// set the print level
	void SetPrintLevel(int n) override;

	// set the solution method
	void SetSolutionMethod(int method);

	// set the zero-initial-guess flag
	void SetZeroInitialGuess(bool b);

private:
	BlockMatrix*			m_pA;		//!< block matrices
	vector<PardisoSolver*>	m_solver;	//!< solvers for solving diagonal blocks

private:
	int		m_method;			//!< 0 = Jacobi, 1 = Gauss-Seidel
	double	m_tol;				//!< convergence tolerance
	int		m_maxiter;			//!< max number of iterations
	int		m_iter;				//!< nr of iterations of last solve
	int		m_printLevel;		//!< set print level
	bool	m_failMaxIter;		//!< fail on max iterations reached
	bool	m_zeroInitGuess;	//!< always use zero as the initial guess

	DECLARE_FECORE_CLASS();
};
