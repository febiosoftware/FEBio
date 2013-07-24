#pragma once

#include "matrix.h"
#include "vector.h"
#include "LinearSolver.h"

class FESolver;

//-----------------------------------------------------------------------------
//! The BFGSSolver solves a nonlinear system of equations using the BFGS method.
//! It depends on the NonLinearSystem to evaluate the function and its jacobian.

class BFGSSolver
{
public:
	//! constructor
	BFGSSolver();

	//! New initialization method
	void Init(int neq, FESolver* pNLS, LinearSolver* pls);

	//! perform a BFGS udpate
	bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1);

	//! solve the equations
	void SolveEquations(vector<double>& x, vector<double>& b);

	//! Performs a linesearch
	double LineSearch(double s);

	//! modified linesearch for Hager-Zhang solver
	double LineSearchCG(double s);

public:
	int		m_maxups;		//!< max nr of QN iters permitted between stiffness reformations
	int		m_maxref;		//!< max nr of reformations per time step
	double	m_cmax;			//!< maximum value for the condition number

	// line search options
	double	m_LSmin;		//!< minimum line search step
	double	m_LStol;		//!< line search tolerance
	int		m_LSiter;		//!< max nr of line search iterations

public:
	// keep a pointer to the linear solver
	LinearSolver*	m_plinsolve;	//!< pointer to linear solver

	// the non-linear system to solve
	FESolver*		m_pNLS;		//!< pointer to nonlinear system to solve
	int				m_neq;		//!< number of equations

	// counters
	int		m_nups;			//!< nr of stiffness updates

	vector<double> m_ui;	//!< displacement increment vector

	// residuals
	vector<double> m_R0;	//!< residual at iteration i-1
	vector<double> m_R1;	//!< residual at iteration i

	// BFGS update vectors
	matrix			m_V;		//!< BFGS update vector
	matrix			m_W;		//!< BFGS update vector
	vector<double>	m_D, m_G, m_H;	//!< temp vectors for calculating BFGS update vectors
};
