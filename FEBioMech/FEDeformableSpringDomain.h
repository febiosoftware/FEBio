#pragma once
#include <FECore/FEDiscreteDomain.h>
#include "FEElasticDomain.h"
#include "FESpringMaterial.h"

//-----------------------------------------------------------------------------
//! domain for deformable springs
class FEDeformableSpringDomain : public FEDiscreteDomain, public FEElasticDomain
{
public:
	//! constructor
	FEDeformableSpringDomain(FEModel* pfem);

	//! Unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat);

	void Activate();

public: // overridden from FEElasticDomain

	//! calculate stiffness matrix
	void StiffnessMatrix(FESolver* psolver);
	void MassMatrix(FESolver* psolver, double scale) {}
	void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) {}

	//! Calculates inertial forces for dynamic problems | todo implement (removed assert DSR)
	void InertialForces(FEGlobalVector& R, vector<double>& F) { }

	//! update domain data
	void Update(const FETimeInfo& tp){}

	//! internal stress forces
	void InternalForces(FEGlobalVector& R);

	//! calculate bodyforces (not used since springs are considered mass-less)
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) {}

protected:
	double InitialLength();
	double CurrentLength();

protected:
	FESpringMaterial*	m_pMat;
	double				m_kbend;	// bending stiffness
	double				m_kstab;	// stabilization penalty
	double				m_L0;	//!< initial spring length

	DECLARE_PARAMETER_LIST();
};
