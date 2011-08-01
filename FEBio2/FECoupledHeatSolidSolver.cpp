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
	//       type of coupled problem. Instead, we'll need to sub-
	//		 matrices which will be stored in the m_Heat and m_Solid
	//		 classes.
//	FESolver::Init();

	// Initialize heat solver
	if (m_Heat.Init() == false) return false;

	// Initialize solid solver
	if (m_Solid.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! Initialize equation system.
//! This solver doesn't manage a linear system, but the two child problems
//! do so we just call the corresponding function.
bool FECoupledHeatSolidSolver::InitEquations()
{
	// Initialize equations for heat problem
	if (m_Heat.InitEquations() == false) return false;

	// Initialize equations for solid problem
	if (m_Solid.InitEquations() == false) return false;

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
void FECoupledHeatSolidSolver::Update(std::vector<double> &u)
{
	// Nothing to do here.
}

//-----------------------------------------------------------------------------
void FECoupledHeatSolidSolver::Serialize(DumpFile &ar)
{
	m_Heat.Serialize(ar);
	m_Solid.Serialize(ar);
}
