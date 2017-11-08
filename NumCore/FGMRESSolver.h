#pragma once
#include "FECore/LinearSolver.h"
#include "CompactMatrix.h"

//-----------------------------------------------------------------------------
//! This class implements an interface to the MKL FGMRES iterative solver for 
//! nonsymmetric indefinite matrices (without pre-conditioning).
class FGMRESSolver : public LinearSolver
{
public:
	//! constructor
	FGMRESSolver();

	//! do any pre-processing (allocates temp storage)
	bool PreProcess();

	//! Factor the matrix (does nothing for iterative solvers)
	bool Factor() { return true; }

	//! Calculate the solution of RHS b and store solution in x
	bool BackSolve(vector<double>& x, vector<double>& b);

	//! Clean up
	void Destroy() {}

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);

	//! Set max nr of iterations
	void SetMaxIterations(int n);

	// Set the print level
	void SetPrintLevel(int n);

	// set residual stopping test flag
	void DoResidualStoppingTest(bool b);

	// set the convergence tolerance for the residual stopping test
	void SetResidualTolerance(double tol);

private:
	int		m_maxiter;			// max nr of iterations
	int		m_print_level;		// output level
	bool	m_doResidualTest;	// do the residual stopping test
	double	m_tol;				// relative residual convergence tolerance

private:
	CompactUnSymmMatrix*	m_pA;		//!< the sparse matrix format
	vector<double>	m_tmp;
};

//-----------------------------------------------------------------------------
//! This class implements an interface to the MKL FGMRES iterative solver with
//! ILUT pre-conditioner for nonsymmetric indefinite matrices.
class FGMRES_ILUT_Solver : public LinearSolver
{
public:
	//! constructor
	FGMRES_ILUT_Solver();

	//! do any pre-processing
	bool PreProcess();

	//! Factor the matrix (does nothing for iterative solvers)
	bool Factor() { return true; }

	//! Calculate the solution of RHS b and store solution in x
	bool BackSolve(vector<double>& x, vector<double>& b);

	//! Clean up
	void Destroy() {}

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);

	//! Set max nr of iterations
	void SetMaxIterations(int n);

	// Set the print level
	void SetPrintLevel(int n);

	// set residual stopping test flag
	void DoResidualStoppingTest(bool b);

	// set the convergence tolerance for the residual stopping test
	void SetResidualTolerance(double tol);

	// set the max fill value
	void SetMaxFill(int n);

	// Set the fill tolerance
	void SetFillTolerance(double fillTol);

private:
	int		m_maxiter;			// max nr of iterations
	int		m_print_level;		// output level
	bool	m_doResidualTest;	// do the residual stopping test
	double	m_tol;				// relative residual convergence tolerance
	int		m_maxfill;
	double	m_fillTol;

private:
	CompactUnSymmMatrix*	m_pA;		//!< the sparse matrix format
	vector<double>	m_tmp;
};

//-----------------------------------------------------------------------------
//! This class implements an interface to the MKL FGMRES iterative solver with
//! ILU0 pre-conditioner for nonsymmetric indefinite matrices.
class FGMRES_ILU0_Solver : public LinearSolver
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
	void Destroy() {}

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);

	//! Set max nr of iterations
	void SetMaxIterations(int n);

	// Set the print level
	void SetPrintLevel(int n);

	// set residual stopping test flag
	void DoResidualStoppingTest(bool b);

	// set the convergence tolerance for the residual stopping test
	void SetResidualTolerance(double tol);

private:
	int		m_maxiter;			// max nr of iterations
	int		m_print_level;		// output level
	bool	m_doResidualTest;	// do the residual stopping test
	double	m_tol;				// relative residual convergence tolerance

private:
	CompactUnSymmMatrix*	m_pA;		//!< the sparse matrix format
	vector<double>	m_tmp;
};
