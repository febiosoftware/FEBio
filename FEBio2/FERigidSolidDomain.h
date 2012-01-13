#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! domain class for 3D rigid elements
//!
class FERigidSolidDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FERigidSolidDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_RIGID_SOLID_DOMAIN; }

	FEDomain* Clone()
	{
		FERigidSolidDomain* pd = new FERigidSolidDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	// update stresses
	void UpdateStresses(FEModel& fem);
};
