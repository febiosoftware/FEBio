#pragma once
#include "FECore/FESolidDomain.h"
#include "FEBioMech/FEElasticDomain.h"
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
//! Domain class for biphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since the biphasic domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBiphasicSolidDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FEBiphasicSolidDomain(FEMesh* pm, FEMaterial* pmat);

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
	void StiffnessMatrix(FESolver* psolver, bool bsymm, double dt);

	//! calculates the global stiffness matrix (steady-state case)
	void StiffnessMatrixSS(FESolver* psolver, bool bsymm, double dt);
	
public:
	// internal work (overridden from FEElasticDomain)
	void InternalForces(FEGlobalVector& R);

	//! internal fluid work
	void InternalFluidWork(vector<double>& R, double dt);

	//! internal fluid work (steady state analysis)
	void InternalFluidWorkSS(vector<double>& R, double dt);

public:
	//! element internal force vector
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);
	
	//! Calculates the internal fluid forces
	bool ElementInternalFluidWork(FESolidElement& elem, vector<double>& fe, double dt);
	
	//! Calculates the internal fluid forces for steady-state response
	bool ElementInternalFluidWorkSS(FESolidElement& elem, vector<double>& fe, double dt);
	
	//! calculates the element biphasic stiffness matrix
	bool ElementBiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm, double dt);
	
	//! calculates the element biphasic stiffness matrix for steady-state response
	bool ElementBiphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm, double dt);
	
	//! calculates the solid element stiffness matrix
	void SolidElementStiffness(FESolidElement& el, matrix& ke);

	//! geometrical stiffness
	void ElementGeometricalStiffness(FESolidElement &el, matrix &ke);

	//! material stiffness component
	void ElementBiphasicMaterialStiffness(FESolidElement& el, matrix& ke);

public: // overridden from FEElasticDomain, but not all implemented in this domain
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf);
    void ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe);
	void InertialForces(FEGlobalVector& R, vector<double>& F) {}
	void StiffnessMatrix(FESolver* psolver) {}
    void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf);
    void ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke);
	void MassMatrix(FESolver* psolver, double scale) {}

protected:
	FEBiphasic*	m_pMat;
};
