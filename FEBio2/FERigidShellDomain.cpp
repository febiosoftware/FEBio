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
//!
void FERigidShellDomain::StiffnessMatrix(FENLSolver* psolver)
{
	FEM& fem = dynamic_cast<FEM&>(psolver->GetFEModel());

	// for dynamic analyses we do need to add the inertial stiffness of the rigid body
	if (fem.GetCurrentStep()->m_nanalysis != FE_DYNAMIC) return;

	matrix ke;

	int NS = m_Elem.size();
	for (int iel=0; iel<NS; ++iel)
	{
		FEShellElement& el = m_Elem[iel];
		assert(el.IsRigid());

		// TODO: add inertial stiffness for shells
	}
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
