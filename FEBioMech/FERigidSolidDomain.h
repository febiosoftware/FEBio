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

	//! Initialize
	bool Initialize(FEModel& fem);

	//! reset data
	void Reset();

public:

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver);

	//! calculates the residual (nothing to do)
	void InternalForces(FESolver* psolver, vector<double>& R);

	//! intertial forces for dynamic problems (overridden from FEElasticDomain)
	void InertialForces(FEGlobalVector& R, vector<double>& F);

	//! body forces (overridden from FEElasticDomain)
	void BodyForce(FEGlobalVector& R, FEBodyForce& BF);

	// update stresses (nothing to do)
	void UpdateStresses(FEModel& fem);

	//! Calculate mass matrix (overridden from FEElasticDomain)
	void MassMatrix(FESolver* psolver, double scale);

	//! body force stiffness (overridden from FEElasticDomain)
	void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf);
};
