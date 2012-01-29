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

	//! Calculates inertial forces for dynamic problems
	void InertialForces(FENLSolver* psolver, vector<double>& R, vector<double>& F) { assert(false); }

	//! calculate residual
//	void Residual(FENLSolver* psolver, vector<double>& R);

	//! update stresses (not used for discrete springs)
	void UpdateStresses(FEModel& fem){}	

	//! internal stress forces
	void InternalForces(FENLSolver* psolver, vector<double>& R);

	//! calculate bodyforces (not used since springs are considered mass-less)
	void BodyForce(FENLSolver* psolver, FEBodyForce& bf, vector<double>& R) {}
};
