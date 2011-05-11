#include "stdafx.h"
#include "FECore/FEDomain.h"
#include "FESolidSolver.h"
#include "FEElasticShellDomain.h"
#include "FEElasticSolidDomain.h"

//=============================================================================
// Rigid Solid
//=============================================================================

//-----------------------------------------------------------------------------
//! Calculates the stiffness matrix for 3D rigid elements
//! Most stiffness contributions are ignored. Only in dynamic problems
//! the inertial stiffness is computed.
//! 
void FERigidSolidDomain::StiffnessMatrix(FESolidSolver* psolver)
{
	FEM& fem = psolver->m_fem;

	// for dynamic analyses we do need to add the inertial stiffness of the rigid body
	if (fem.m_pStep->m_nanalysis != FE_DYNAMIC) return;

	// element stiffness matrix
	matrix ke;

	// repeat over all solid elements
	int NE = m_Elem.size();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];
		assert(el.IsRigid());

		UnpackElement(el);

		int ndof = 3*el.Nodes();
		ke.Create(ndof, ndof);
		ke.zero();

		// add the inertial stiffness for dynamics
		ElementInertialStiffness(fem, el, ke);

		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, el.LM(), ke);
	}
}

//-----------------------------------------------------------------------------
//! Calculates the residual vector 3D rigid elements
//! we only calculate the body-forces for rigid elements
//!
void FERigidSolidDomain::Residual(FESolidSolver *psolver, vector<double>& R)
{
	FEM& fem = psolver->m_fem;

	if (fem.HasBodyForces() == false) return;

	// element force vector
	vector<double> fe;

	// loop over all elements
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];
		assert(el.IsRigid());

		// unpack the element
		UnpackElement(el);

		FEMaterial* pm = fem.GetMaterial(el.GetMatID());

		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);

		// apply body force to rigid elements
		BodyForces(fem, el, fe);

		// assemble element 'fe'-vector into global R vector
		psolver->AssembleResidual(el.m_node, el.LM(), fe, R);
	}
}

//-----------------------------------------------------------------------------
//! We don't need to update the stresses for rigid elements
//!
void FERigidSolidDomain::UpdateStresses(FEM &fem)
{
	// Nothing to see here. Please move on.
}

//=============================================================================
// Rigid Shell
//=============================================================================


//-----------------------------------------------------------------------------
//! Calculate stiffness contributions for rigid shells.
//!
void FERigidShellDomain::StiffnessMatrix(FESolidSolver* psolver)
{
	FEM& fem = psolver->m_fem;

	// for dynamic analyses we do need to add the inertial stiffness of the rigid body
	if (fem.m_pStep->m_nanalysis != FE_DYNAMIC) return;

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
void FERigidShellDomain::Residual(FESolidSolver* psolver, vector<double>& R)
{
	FEM& fem = psolver->m_fem;

	if (fem.HasBodyForces() == false) return;

	// element force vector
	vector<double> fe;

	int NS = m_Elem.size();
	for (int i=0; i<NS; ++i)
	{
		// get the element
		FEShellElement& el = m_Elem[i];

		assert(el.IsRigid());

		// create the element force vector and initialize to zero
		int ndof = 6*el.Nodes();
		fe.assign(ndof, 0);

		UnpackElement(el);

		// apply body forces to shells
		BodyForces(fem, el, fe);

		// assemble the residual
		psolver->AssembleResidual(el.m_node, el.LM(), fe, R);
	}
}

//-----------------------------------------------------------------------------
//! update stresses for rigid shells.
//!
void FERigidShellDomain::UpdateStresses(FEM &fem)
{
	// Nothing to see here. Please move on.
}
