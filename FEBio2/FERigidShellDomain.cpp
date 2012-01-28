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
void FERigidShellDomain::Residual(FENLSolver* psolver, vector<double>& R)
{
	FEM& fem = dynamic_cast<FEM&>(psolver->GetFEModel());

	if (fem.HasBodyForces() == false) return;

	// element force vector
	vector<double> fe;

	vector<int> lm;

	int NS = m_Elem.size();
	for (int i=0; i<NS; ++i)
	{
		// get the element
		FEShellElement& el = m_Elem[i];

		assert(el.IsRigid());

		// create the element force vector and initialize to zero
		int ndof = 6*el.Nodes();
		fe.assign(ndof, 0);

		// apply body forces to shells
		ElementBodyForce(fem, el, fe);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble the residual
		psolver->AssembleResidual(el.m_node, lm, fe, R);
	}
}

//-----------------------------------------------------------------------------
//! update stresses for rigid shells.
//!
void FERigidShellDomain::UpdateStresses(FEModel &fem)
{
	// Nothing to see here. Please move on.
}
