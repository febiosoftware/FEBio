#include "stdafx.h"
#include "FESolidSolver.h"
#include "FEElasticSolidDomain.h"

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolidSolver::Augment
//  This functions performs the Lagrange augmentations
//  It returns true if all the augmentation have converged, 
//	otherwise it returns false
//
// TODO: There is an inherent problem with this approach. Since
//  Lagrangian multipliers are inherited from previous timesteps
//  they might not be zero in case a node-surface contact breaks. 
//  The node's gap value needs to become negative to a certain value
//  before the Lagr. multipliers dissapears. 
//

bool FESolidSolver::Augment()
{
	// Assume we will pass (can't hurt to be optimistic)
	bool bconv = true;

	// Do rigid joint augmentations
	if (!m_fem.m_RJ.empty())
	{
		// loop over all rigid joints
		for (int i=0; i<(int) m_fem.m_RJ.size(); ++i) bconv = m_fem.m_RJ[i]->Augment() && bconv;
	}

	// Do contact augmentations
	if (m_fem.ContactInterfaces() > 0)
	{
		// loop over all contact interfaces
		for (int i=0; i<m_fem.ContactInterfaces(); ++i) bconv = m_fem.ContactInterface(i)->Augment(m_naug) && bconv;
	}

	// do linear constraint augmentations
	if (m_fem.m_LCSet.size())
	{
		int n = m_fem.m_LCSet.size();
		list<FELinearConstraintSet*>::iterator im = m_fem.m_LCSet.begin();
		for (int i=0; i<n; ++i, ++im) bconv = (*im)->Augment(m_naug) && bconv;
	}

	// do incompressibility multipliers for 3Field domains
	FEMesh& mesh = m_fem.GetMesh();
	int ND = mesh.Domains();
	for (int i=0; i<ND; ++i)
	{
		FE3FieldElasticSolidDomain* pd = dynamic_cast<FE3FieldElasticSolidDomain*>(&mesh.Domain(i));
		if (pd) bconv = (pd->Augment() && bconv);
	}

	return bconv;
}
