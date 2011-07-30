#include "stdafx.h"
#include "FECoupledHeatSolidSolver.h"

//-----------------------------------------------------------------------------
//! class constructor
FECoupledHeatSolidSolver::FECoupledHeatSolidSolver(FEM& fem) : FESolver(fem), m_Heat(fem), m_Solid(fem)
{
	m_niter = 0;
}

//-----------------------------------------------------------------------------
//! Initialization
bool FECoupledHeatSolidSolver::Init()
{
	// NOTE: Don't call base class since it will try to allocate
	//       a global stiffness matrix. We don't need one for this
	//       type of coupled problem.
//	FESolver::Init();

	// TODO: The solvers use FEM::m_neq to determine the size
	//       for the solution vectors. Obviously that won't 
	//       work here. I have to figure out a different way
	//       to determine equation numbers.

	// Initialize heat solver
	if (m_Heat.Init() == false) return false;

	// Initialize solid solver
	if (m_Solid.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
bool FECoupledHeatSolidSolver::SolveStep(double time)
{
	// First we solve the heat problem
	if (m_Heat.SolveStep(time) == false) return false;

	// Now we project the nodal temperatures
	// to the integration points, and use them to set up an
	// initial stress for the linear solid solver
	// TODO: do above

	// Next, we solve the linear solid problem
	if (m_Solid.SolveStep(time) == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FECoupledHeatSolidSolver::Serialize(DumpFile &ar)
{
	m_Heat.Serialize(ar);
	m_Solid.Serialize(ar);
}
