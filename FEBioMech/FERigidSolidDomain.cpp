#include "stdafx.h"
#include "FERigidSolidDomain.h"

//-----------------------------------------------------------------------------
FEDomain* FERigidSolidDomain::Clone()
{
	FERigidSolidDomain* pd = new FERigidSolidDomain(m_pMesh, m_pMat);
	pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
	return pd;
}

//-----------------------------------------------------------------------------
// We need to override it since the base class version will not work for rigid domains.
void FERigidSolidDomain::Reset()
{
	// nothing to reset here
}

//-----------------------------------------------------------------------------
//! Calculates the stiffness matrix for 3D rigid elements.
//! Rigid elements don't generate stress, so there is nothing to do here
void FERigidSolidDomain::StiffnessMatrix(FESolver* psolver)
{
	// Caught you looking!
}

//-----------------------------------------------------------------------------
// Rigid bodies do not generate stress so there is nothing to do here
void FERigidSolidDomain::InternalForces(FESolver* psolver, vector<double>& R)
{
	// what you looking at ?!
}

//-----------------------------------------------------------------------------
//! We don't need to update the stresses for rigid elements
//!
void FERigidSolidDomain::UpdateStresses(FEModel &fem)
{
	// Nothing to see here. Please move on.
}
