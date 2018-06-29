#pragma once
#include <FECore/LinearSolver.h>
#include "CompactUnSymmMatrix.h"
#include "Preconditioner.h"

//-----------------------------------------------------------------------------
//! This class implements an interface to the MKL FGMRES iterative solver with
//! ILU0 pre-conditioner for nonsymmetric indefinite matrices.
class FGMRES_ILU0_Solver : public IterativeLinearSolver
{
public:
	//! constructor
	FGMRES_ILU0_Solver();

	//! do any pre-processing
	bool PreProcess();

	//! Factor the matrix (does nothing for iterative solvers)
	bool Factor() { return true; }

	//! Calculate the solution of RHS b and store solution in x
	bool BackSolve(vector<double>& x, vector<double>& b);

	//! Clean up
	void Destroy() override;

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);

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

	// set the convergence tolerance for the residual stopping test
	void SetResidualTolerance(double tol);

	// do the zero diagonal check during preconditioner
	void DoZeroDiagonalCheck(bool b);

	// Set the zero diagonal tolerance value
	void SetZeroDiagonalTolerance(double tol);

	// set the zero diagonal replacement value
	void SetZeroDiagonalReplacement(double val);

	//! This solver uses a preconditioner
	bool HasPreconditioner() const override { return true; }

private:
	int		m_maxiter;			// max nr of iterations
	int		m_nrestart;			// max nr of non-restarted iterations
	int		m_print_level;		// output level
	bool	m_doResidualTest;	// do the residual stopping test
	double	m_tol;				// relative residual convergence tolerance

private:
	SparseMatrix*			m_pA;		//!< the sparse matrix format
	ILU0_Preconditioner		m_P;		//!< the preconditioner
	vector<double>	m_tmp;
	bool			m_doPreCond;
};
