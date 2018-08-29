#pragma once

#include "matrix.h"
#include "vector.h"
#include "LinearSolver.h"
#include "FENewtonStrategy.h"

//-----------------------------------------------------------------------------
//! The BFGSSolver solves a nonlinear system of equations using the BFGS method.
//! It depends on the NonLinearSystem to evaluate the function and its jacobian.

class FECORE_API BFGSSolver : public FENewtonStrategy
{
public:
	//! constructor
	BFGSSolver(FENewtonSolver* pns);

	//! New initialization method
	void Init(int neq, LinearSolver* pls);

	//! perform a BFGS udpate
	bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1);

	//! solve the equations
	void SolveEquations(vector<double>& x, vector<double>& b);

public:
	// keep a pointer to the linear solver
	LinearSolver*	m_plinsolve;	//!< pointer to linear solver
	int				m_neq;		//!< number of equations

	// BFGS update vectors
	matrix			m_V;		//!< BFGS update vector
	matrix			m_W;		//!< BFGS update vector
	vector<double>	m_D, m_G, m_H;	//!< temp vectors for calculating BFGS update vectors
};
