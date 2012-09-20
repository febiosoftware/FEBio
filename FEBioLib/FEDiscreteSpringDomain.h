#pragma once
#include "FECore/FEDiscreteDomain.h"
#include "FEElasticDomain.h"

//-----------------------------------------------------------------------------
//! domain for discrete elements
class FEDiscreteSpringDomain : public FEDiscreteDomain, public FEElasticDomain
{
public:
	//! constructor
	FEDiscreteSpringDomain(FEMesh* pm, FEMaterial* pmat) : FEDiscreteDomain(FE_DISCRETE_DOMAIN, pm, pmat) {}

	//! Clone this domain
	FEDomain* Clone();

	//! Unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! Serialize data to archive
	void Serialize(DumpFile& ar);

public: // overridden from FEElasticDomain

	//! calculate stiffness matrix
	void StiffnessMatrix(FENLSolver* psolver);
	void InertialStiffness   (FENLSolver* psolver) {}
	void BodyForceStiffness  (FENLSolver* psolver, FEBodyForce& bf) {}

	//! Calculates inertial forces for dynamic problems
	void InertialForces(FEGlobalVector& R, vector<double>& F) { assert(false); }

	//! calculate residual
//	void Residual(FEGlobalVector& R);

	//! update stresses (not used for discrete springs)
	void UpdateStresses(FEModel& fem){}	

	//! internal stress forces
	void InternalForces(FEGlobalVector& R);

	//! calculate bodyforces (not used since springs are considered mass-less)
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) {}
};
