#include "stdafx.h"
#include "FERigidShellDomain.h"
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
FEDomain* FERigidShellDomain::Clone()
{
	FERigidShellDomain* pd = new FERigidShellDomain(m_pMesh, m_pMat);
	pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
	return pd;
}

//-----------------------------------------------------------------------------
//! Calculate stiffness contributions for rigid shells.
//! Since rigid elements don't generate stress, we don't need to do
//! anything here.
void FERigidShellDomain::StiffnessMatrix(FENLSolver* psolver)
{
	// Caught you looking!
}


//-----------------------------------------------------------------------------
//! calculate residual forces for rigid shells
//!
void FERigidShellDomain::InternalForces(FENLSolver* psolver, vector<double>& R)
{
	// Nothing to do.
}

//-----------------------------------------------------------------------------
//! update stresses for rigid shells.
//!
void FERigidShellDomain::UpdateStresses(FEModel &fem)
{
	// Nothing to see here. Please move on.
}
