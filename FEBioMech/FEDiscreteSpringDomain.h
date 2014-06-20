#pragma once
#include "FECore/FEDiscreteDomain.h"
#include "FEElasticDomain.h"
#include "FESpringMaterial.h"

//-----------------------------------------------------------------------------
//! domain for discrete elements
class FEDiscreteSpringDomain : public FEDiscreteDomain, public FEElasticDomain
{
public:
	//! constructor
	FEDiscreteSpringDomain(FEMesh* pm, FEMaterial* pmat);

	//! Unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

public: // overridden from FEElasticDomain

	//! calculate stiffness matrix
	void StiffnessMatrix(FESolver* psolver);
	void MassMatrix(FESolver* psolver, double scale) {}
	void BodyForceStiffness  (FESolver* psolver, FEBodyForce& bf) {}

	//! Calculates inertial forces for dynamic problems
	void InertialForces(FESolver* psolver, FEGlobalVector& R, vector<double>& F) { assert(false); }

	//! update stresses (not used for discrete springs)
	void UpdateStresses(FEModel& fem){}	

	//! internal stress forces
	void InternalForces(FEGlobalVector& R);

	//! calculate bodyforces (not used since springs are considered mass-less)
	void BodyForce(FESolver* psolver, FEGlobalVector& R, FEBodyForce& bf) {}

protected:
	FESpringMaterial*	m_pMat;
};
