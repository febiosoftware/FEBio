#pragma once
#include <FECore/LinearSolver.h>
#include <FECore/SparseMatrix.h>
#include <FECore/Preconditioner.h>

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

	//! Factor the matrix (does nothing for iterative solvers)
	bool Factor() override { return true; }

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

	//! Set the nr of non-restarted iterations
	void SetNonRestartedIterations(int n);

	// Set the print level
	void SetPrintLevel(int n);

	// set residual stopping test flag
	void DoResidualStoppingTest(bool b);

	// set zero norm stopping test flag
	void DoZeroNormStoppingTest(bool b);

	// set the convergence tolerance for the residual stopping test
	void SetResidualTolerance(double tol);

	//! This solver does not use a preconditioner
	bool HasPreconditioner() const override;

	//! convenience function for solving linear system Ax = b
	bool Solve(SparseMatrix* A, vector<double>& x, vector<double>& b);

public:
	// set the preconditioner
	void SetPreconditioner(Preconditioner* P) override;

private:
	int		m_maxiter;			// max nr of iterations
	int		m_nrestart;			// max nr of non-restarted iterations
	int		m_print_level;		// output level
	bool	m_doResidualTest;	// do the residual stopping test
	bool	m_doZeroNormTest;	// do the zero-norm stopping test
	double	m_tol;				// relative residual convergence tolerance

private:
	SparseMatrix*	m_pA;		//!< the sparse matrix format
	Preconditioner*	m_P;		//!< the preconditioner
	vector<double>	m_tmp;
	bool			m_doPreCond;
};
