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
#include <FECore/SparseMatrix.h>

//-----------------------------------------------------------------------------
//! This class implements an interface to the MKL FGMRES iterative solver for 
//! nonsymmetric indefinite matrices (without pre-conditioning).
class FGMRESSolver : public IterativeLinearSolver
{
public:
	//! constructor
	FGMRESSolver(FEModel* fem);

	//! do any pre-processing (allocates temp storage)
	bool PreProcess() override;

	//! Factor the matrix
	bool Factor() override;

	//! Calculate the solution of RHS b and store solution in x
	bool BackSolve(double* x, double* b) override;

	//! Clean up
	void Destroy() override;

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	//! Set the sparse matrix
	bool SetSparseMatrix(SparseMatrix* pA) override;

	//! Set max nr of iterations
	void SetMaxIterations(int n);

	//! Get the max nr of iterations
	int GetMaxIterations() const;

	//! Set the nr of non-restarted iterations
	void SetNonRestartedIterations(int n);

	// Set the print level
	void SetPrintLevel(int n) override;

	// set residual stopping test flag
	void DoResidualStoppingTest(bool b);

	// set zero norm stopping test flag
	void DoZeroNormStoppingTest(bool b);

	// set the relative convergence tolerance for the residual stopping test
	void SetRelativeResidualTolerance(double tol);

	// set the absolute convergence tolerance for the residual stopping test
	void SetAbsoluteResidualTolerance(double tol);

	//! This solver does not use a preconditioner
	bool HasPreconditioner() const override;

	//! convenience function for solving linear system Ax = b
	bool Solve(SparseMatrix* A, vector<double>& x, vector<double>& b);

	//! fail if max iterations reached
	void FailOnMaxIterations(bool b);

	//! print the condition number
	void PrintConditionNumber(bool b);

	// do jacobi preconditioning
	void DoJacobiPreconditioning(bool b);

public:
	// set the preconditioner
	void SetLeftPreconditioner(LinearSolver* P) override;
	void SetRightPreconditioner(LinearSolver* P) override;

	// get the preconditioner
	LinearSolver* GetLeftPreconditioner() override;
	LinearSolver* GetRightPreconditioner() override;

protected:
	SparseMatrix* GetSparseMatrix() { return m_pA; }

private:
	int		m_maxiter;			// max nr of iterations
	int		m_nrestart;			// max nr of non-restarted iterations
	int		m_print_level;		// output level
	bool	m_doResidualTest;	// do the residual stopping test
	bool	m_doZeroNormTest;	// do the zero-norm stopping test
	double	m_reltol;			// relative residual convergence tolerance
	double	m_abstol;			// absolute residual tolerance
	bool	m_maxIterFail;
	bool	m_print_cn;			// Calculate and print the condition number
	bool	m_do_jacobi;

private:
	SparseMatrix*	m_pA;		//!< the sparse matrix format
	LinearSolver*	m_P;		//!< the left preconditioner
	LinearSolver*	m_R;		//!< the right preconditioner
	vector<double>	m_tmp;
	vector<double>	m_Rv;		//!< used when a right preconditioner is ued
	vector<double>	m_W;		//!< Jacobi preconditioner

	DECLARE_FECORE_CLASS();
};
