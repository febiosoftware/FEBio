#pragma once
#include "FECore/FESolidDomain.h"
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"
#include <FECore/FETypes.h>

//-----------------------------------------------------------------------------
//! domain described by Lagrange-type 3D volumetric elements
//!
class FEElasticSolidDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FEElasticSolidDomain(FEModel* pfem);

	//! assignment operator
	FEElasticSolidDomain& operator = (FEElasticSolidDomain& d);

	//! initialize class
	bool Initialize(FEModel& fem);

	//! activate
	void Activate();

	//! initialize elements
	virtual void InitElements();

	//! Unpack solid element data
	void UnpackLM(FEElement& el, vector<int>& lm);

public: // overrides from FEDomain

	//! get the material
	FEMaterial* GetMaterial() { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pm);

public: // overrides from FEElasticDomain

	// update stresses
	void Update(const FETimePoint& tp);

	// update the element stress
	virtual void UpdateElementStress(int iel, double dt);

	//! intertial forces for dynamic problems
	void InertialForces(FEGlobalVector& R, vector<double>& F);

	//! intertial forces for dynamic problems (used by FESolidSolver2)
	void InertialForces2(FEGlobalVector& R, vector<double>& F);

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
	virtual void ElementStiffness(const FETimePoint& tp, int iel, matrix& ke);

	//! geometrical stiffness (i.e. initial stress)
	virtual void ElementGeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness component
	virtual void ElementMaterialStiffness(FESolidElement& el, matrix& ke);

	//! calculates the solid element mass matrix
	void ElementMassMatrix(FESolidElement& el, matrix& ke, double a);

	//! calculates the stiffness matrix due to body forces 
	void ElementBodyForceStiffness(FEBodyForce& bf, FESolidElement& el, matrix& ke);

	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for solid elements
	virtual void ElementInternalForce(FESolidElement& el, vector<double>& fe);

	//! Calculatess external body forces for solid elements
	void ElementBodyForce(FEBodyForce& BF, FESolidElement& elem, vector<double>& fe);

	// ---

protected:
	FESolidMaterial*	m_pMat;
};
