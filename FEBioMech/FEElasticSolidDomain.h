#pragma once
#include "FECore/FESolidDomain.h"
#include "FECore/FEModel.h"
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! domain described by Lagrange-type 3D volumetric elements
//!
class FEElasticSolidDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FEElasticSolidDomain(FEMesh* pm, FEMaterial* pmat);

	//! \todo Do I really use this?
	FEElasticSolidDomain& operator = (FEElasticSolidDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! initialize class
	bool Initialize(FEModel& fem);

	//! initialize elements
	void InitElements();

	//! Unpack solid element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

public: // overrides from FEElasticDomain

	// update stresses
	void UpdateStresses(FEModel& fem);

	// update the element stress
	void UpdateElementStress(int iel, double dt);

	//! intertial forces for dynamic problems
	void InertialForces(FEGlobalVector& R, vector<double>& F);

	//! internal stress forces
	void InternalForces(FEGlobalVector& R);

	//! body forces
	void BodyForce(FEGlobalVector& R, FEBodyForce& BF);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver);

	//! calculates inertial stiffness
	void MassMatrix(FESolver* psolver, double scale);

	//! body force stiffness
	void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf);

public:
	// --- S T I F F N E S S ---

	//! calculates the solid element stiffness matrix
	virtual void ElementStiffness(FEModel& fem, int iel, matrix& ke);

	//! geometrical stiffness (i.e. initial stress)
	void ElementGeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness component
	void ElementMaterialStiffness(FESolidElement& el, matrix& ke);

	//! calculates the solid element mass matrix
	void ElementMassMatrix(FESolidElement& el, matrix& ke, double a);

	//! calculates the stiffness matrix due to body forces 
	void ElementBodyForceStiffness(FEBodyForce& bf, FESolidElement& el, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

	//! Calculatess external body forces for solid elements
	void ElementBodyForce(FEBodyForce& BF, FESolidElement& elem, vector<double>& fe);

	// ---

protected:
	FESolidMaterial*	m_pMat;
};
