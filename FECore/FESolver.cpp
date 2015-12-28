#include "stdafx.h"
#include "FESolver.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
FESolver::FESolver(FEModel* pfem) : FECoreBase(FESOLVER_ID), m_fem(*pfem)
{ 
	m_bsymm = true; // assume symmetric stiffness matrix
	m_niter = 0;
}

//-----------------------------------------------------------------------------
FESolver::~FESolver()
{
}

//-----------------------------------------------------------------------------
//! Get the FE model
FEModel& FESolver::GetFEModel()
{ 
	return m_fem; 
}

//-----------------------------------------------------------------------------
//! This function is called right before SolveStep and should be used to initialize
//! time dependent information and other settings.
bool FESolver::InitStep(double time)
{
	FEModel& fem = GetFEModel();

	// evaluate load curve values at current time
	fem.EvaluateLoadCurves(time);

	// evaluate the parameter lists
	fem.EvaluateAllParameterLists();

	// re-initialize materials
	// TODO: I need to do this since the material parameters can have changed and thus a new initialization
	//       needs to be done to see if the material parameters are still valid. I would like to add value checking
	//       directly in the parameter evaluation above so this can be removed.
	if (fem.InitMaterials() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
bool FESolver::Init()
{
	return true;
}
