#include "stdafx.h"
#include "FENLSolver.h"

//-----------------------------------------------------------------------------
FENLSolver::FENLSolver(FEModel& fem) : m_fem(fem)
{ 
	m_bsymm = true; // assume symmetric stiffness matrix
	m_solvertype = 0;
}

//-----------------------------------------------------------------------------
FENLSolver::~FENLSolver()
{
}

//-----------------------------------------------------------------------------
//! Get the FE model
FEModel& FENLSolver::GetFEModel()
{ 
	return m_fem; 
}
