#pragma once

#include "FECore/SparseMatrix.h"
#include "FECore/LinearSolver.h"

//-----------------------------------------------------------------------------
//! The NonLinearSystem class describes a system of nonlinear equations. It
//! provides two functions: eval() which takes a point and calculates
//! the function value at this point; and jacobian() which calculates the jacobian 
//! of the system at a particular point.

class NonLinearSystem
{
public:
	//! constructor
	NonLinearSystem(int neq, int ndof);

	//! destructor
	virtual ~NonLinearSystem();

	//! evaluate the system and its jacobian at the point x
	virtual void eval(vector<double>& x, vector<double>& y, SparseMatrix* pK = 0) = 0;

	//! return the number of nonlinear equations
	int equations() { return m_neq; }

	//! return the number of degrees of freedom
	int dofs() { return m_ndof; }

protected:
	int	m_neq;	//!< nr of equations
	int	m_ndof;	//!< nr of degrees of freedom
};

//-----------------------------------------------------------------------------
//! The BFGSSolver solves a nonlinear system of equations using the BFGS method.
//! It depends on the NonLinearSystem to evaluate the function and its jacobian.

class BFGSSolver
{
public:
	//! constructor
	BFGSSolver();

	//! destructor
	virtual ~BFGSSolver();

	//! solve the nonlinear system of equations
	bool Solve(vector<double>& x, NonLinearSystem& S);

protected:
	//! check convergence
	virtual bool converged();

	//! find the next search direction
	void findsd();

	//! do a line search
	void linesearch();

	//! update the solution vector
	void update();

protected:
	NonLinearSystem*	m_pS;	//!< pointer to the nonlinear system being solved
	SparseMatrix*		m_pK;	//!< pointer to the jacobian of the nonlinear system
	LinearSolver*		m_pls;	//!< pointer to a linear solver

	vector<double>	m_x;	//!< solution vector
	vector<double>	m_xt;	//!< trial solution vector (used in line search)
	vector<double>	m_u;	//!< search direction
	vector<double>	m_F0;	//!< function value at previous solution
	vector<double>	m_F1;	//!< function value at current solution

	vector<double>	m_D;	//!< help vector
	vector<double>	m_G;	//!< help vector
	vector<double>	m_H;	//!< help vector

	vector< vector<double> >	m_V;	//!< BFGS update vector v
	vector< vector<double> >	m_W;	//!< BFGS update vector w

	double	m_ls;	//!< line search factor

	int	m_nups;		//!< nr of BFGS updates
	double	m_cmax;	//!< max condition number

	int	m_maxups;	//!< max nr of updates
};
