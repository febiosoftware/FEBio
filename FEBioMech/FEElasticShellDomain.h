#pragma once
#include "FECore/FEShellDomain.h"
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticShellDomain : public FEShellDomain, public FEElasticDomain
{
public:
	FEElasticShellDomain(FEMesh* pm, FEMaterial* pmat);

	//! \todo do I really need this?
	FEElasticShellDomain& operator = (FEElasticShellDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! Initialize domain
	bool Initialize(FEModel& fem);

	//! Unpack shell element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

public: // overrides from FEElasticDomain

	//! calculates the residual
//	void Residual(FESolver* psolver, vector<double>& R);

	//! internal stress forces
	void InternalForces(FEGlobalVector& R);

	//! Calculates inertial forces for dynamic problems
	void InertialForces(FEGlobalVector& R, vector<double>& F) { assert(false); }

	//! calculate body force
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf);

	// update stresses
	void UpdateStresses(FEModel& fem);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver);

	// inertial stiffness \todo implement this
	void MassMatrix(FESolver* psolver, double scale) { assert(false); }

	// body force stiffness \todo implement this
	void BodyForceStiffness  (FESolver* psolver, FEBodyForce& bf) { assert(false); }

public:

	// --- S T I F F N E S S --- 

	//! calculates the shell element stiffness matrix
	void ElementStiffness(int iel, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for shell elements
	void ElementInternalForce(FEShellElement& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEModel& fem, FEShellElement& el, vector<double>& fe);

	//! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe);

protected:
	FESolidMaterial*	m_pMat;
};
