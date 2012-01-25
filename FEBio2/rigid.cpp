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
//! Calculates the residual vector 3D rigid elements
//! we only calculate the body-forces for rigid elements
//!
void FERigidSolidDomain::Residual(FENLSolver *psolver, vector<double>& R)
{
	FEM& fem = dynamic_cast<FEM&>(psolver->GetFEModel());

	if (fem.HasBodyForces() == false) return;

	// element force vector
	vector<double> fe;

	vector<int> lm;

	// loop over all elements
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
		assert(el.IsRigid());

		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);

		// apply body force to rigid elements
		BodyForces(fem, el, fe);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element 'fe'-vector into global R vector
		psolver->AssembleResidual(el.m_node, lm, fe, R);
	}
}

//-----------------------------------------------------------------------------
//! We don't need to update the stresses for rigid elements
//!
void FERigidSolidDomain::UpdateStresses(FEModel &fem)
{
	// Nothing to see here. Please move on.
}
