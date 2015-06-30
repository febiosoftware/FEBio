#pragma once
#include <FECore\FESolidDomain.h>
#include <FEBioMech\FEElasticDomain.h>
#include "FEThermoElasticMaterial.h"

//-----------------------------------------------------------------------------
// This class implements a thermo-elastic domain.
class FEThermoElasticSolidDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FEThermoElasticSolidDomain(FEMesh* pm, FEMaterial* pmat);

	//! initialize class
	bool Initialize(FEModel& fem);

	//! reset domain data
	void Reset();

	//! intitialize element data
	void InitElements();

	//! Unpack solid element data  (overridden from FEDomain)
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

public:
	// update stresses (overridden from FEElasticDomain)
	void UpdateStresses(FEModel& fem);

	// update element stress
	void UpdateElementStress(int iel);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver);

public:
	// internal work (overridden from FEElasticDomain)
	void InternalForces(FEGlobalVector& R);

	//! internal thermal work
	void InternalThermalWork(vector<double>& R);

public:
	//! element internal force vector
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);
	
	//! Calculates the internal temperature forces
	bool ElementInternalThermalWork(FESolidElement& elem, vector<double>& fe);
	
	//! calculates the element biphasic stiffness matrix
	void ElementStiffness(FESolidElement& el, matrix& ke);
	
	//! calculates the solid element stiffness matrix
	void SolidElementStiffness(FESolidElement& el, matrix& ke);

	//! geometrical stiffness
	void ElementGeometricalStiffness(FESolidElement &el, matrix &ke);

	//! material stiffness
	void ElementMaterialStiffness(FESolidElement &el, matrix &ke);

	//! Conduction stiffness
	void ElementConductionStiffness(FESolidElement &el, matrix& ke);

	//! Thermal stiffness
	void ElementThermalStiffness(FESolidElement &el, matrix& ke);

	//! Conductivity gradient stiffness
	void ElementGradientStiffness(FESolidElement &el, matrix& ke);

public: // overridden from FEElasticDomain, but not all implemented in this domain
	void InertialForces(FEGlobalVector& R, vector<double>& F) {}
	void MassMatrix(FESolver* psolver, double scale) {}
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) {}
	void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) {}

protected:
	FEThermoElasticMaterial*	m_pMat;
};
