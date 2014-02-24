#include "stdafx.h"
#include "FESolver.h"

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
