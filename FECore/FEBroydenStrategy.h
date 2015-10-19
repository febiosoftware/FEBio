#pragma once
#include "matrix.h"
#include "FENewtonStrategy.h"

//-----------------------------------------------------------------------------
//! This class implements the Broyden quasi-newton strategy. 
class FEBroydenStrategy : public FENewtonStrategy
{
public:
	//! constructor
	FEBroydenStrategy();

	//! Initialization
	void Init(int neq, LinearSolver* pls);

	//! perform a quasi-Newton udpate
	bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1);

	//! solve the equations
	void SolveEquations(vector<double>& x, vector<double>& b);

private:
	// keep a pointer to the linear solver
	LinearSolver*	m_plinsolve;	//!< pointer to linear solver
	int				m_neq;			//!< number of equations

	// Broyden update vectors
	matrix			m_R;		//!< Broyden update vector "r"
	matrix			m_D;		//!< Broydeb update vector "delta"
	vector<double>	m_rho;		//!< temp vectors for calculating Broyden update vectors
};
