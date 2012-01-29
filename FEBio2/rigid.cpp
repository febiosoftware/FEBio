#include "stdafx.h"
#include "FECore/FEDomain.h"
#include "FESolidSolver.h"
#include "FERigidSolidDomain.h"

//=============================================================================
// Rigid Solid
//=============================================================================

//-----------------------------------------------------------------------------
//! Calculates the stiffness matrix for 3D rigid elements
//! Most stiffness contributions are ignored. Only in dynamic problems
//! the inertial stiffness is computed.
//! 
void FERigidSolidDomain::StiffnessMatrix(FENLSolver* psolver)
{
	FEM& fem = dynamic_cast<FEM&>(psolver->GetFEModel());

	// for dynamic analyses we do need to add the inertial stiffness of the rigid body
	if (fem.GetCurrentStep()->m_nanalysis != FE_DYNAMIC) return;

	// element stiffness matrix
	matrix ke;

	vector<int> lm;

	// repeat over all solid elements
	int NE = m_Elem.size();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];
		assert(el.IsRigid());

		int ndof = 3*el.Nodes();
		ke.resize(ndof, ndof);
		ke.zero();

		// add the inertial stiffness for dynamics
		ElementInertialStiffness(fem, el, ke);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
// Rigid bodies do not generate stress so there is nothing to do here
void FERigidSolidDomain::InternalForces(FENLSolver* psolver, vector<double>& R)
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
