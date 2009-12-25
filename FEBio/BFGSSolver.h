#pragma once

#include "FECore/matrix.h"
#include "FECore/vector.h"
#include "FECore/LinearSolver.h"

//-----------------------------------------------------------------------------
//! The BFGSSolver solves a nonlinear system of equations using the BFGS method.
//! It depends on the NonLinearSystem to evaluate the function and its jacobian.

class BFGSSolver
{
public:
	//! constructor
	BFGSSolver();

	//! initialization
	void Init(int neq, LinearSolver* pls);

	//! perform a BFGS udpate
	bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1);

	//! solve the equations
	void SolveEquations(vector<double>& x, vector<double>& b);

public:
	int		m_maxups;		//!< max nr of QN iters permitted between stiffness reformations
	int		m_maxref;		//!< max nr of reformations per time step

	// keep a pointer to the linear solver
	LinearSolver*	m_plinsolve;	//!< pointer to linear solver

	// counters
	int		m_nups;			//!< nr of stiffness updates
	double	m_cmax;			//!< maximum value for the condition number

	// BFGS update vectors
	matrix			m_V;
	matrix			m_W;
	vector<double>	m_D, m_G, m_H;	//!< temp vectors for calculating BFGS update vectors
};
