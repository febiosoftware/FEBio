#pragma once
#include "FEElasticShellDomain.h"

//-----------------------------------------------------------------------------
//! domain class for 3D rigid shell elements
//!
class FERigidShellDomain : public FEElasticShellDomain
{
public:
	//! constructor
	FERigidShellDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticShellDomain(pm, pmat) { m_ntype = FE_RIGID_SHELL_DOMAIN; }

	//! clone this domain
	FEDomain* Clone();

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolidSolver* psolver);

	//! calculates the residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	// update stresses
	void UpdateStresses(FEModel& fem);
};
