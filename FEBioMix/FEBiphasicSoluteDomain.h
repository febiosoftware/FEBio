#pragma once
#include "FECore/FESolidDomain.h"
#include "FEBiphasicSolute.h"
#include "FEBioMech/FEElasticDomain.h"
#include "FECore/FETypes.h"

//-----------------------------------------------------------------------------
//! Domain class for biphasic-solute 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since this domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBiphasicSoluteDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FEBiphasicSoluteDomain(FEModel* pfem);
	
	//! reset domain data
	void Reset();

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat);

	//! Unpack solid element data (overridden from FEDomain)
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! initialize class
	bool Initialize(FEModel& fem);

	//! Activate
	void Activate();

	//! initialize elements for this domain
	void InitElements();

	// update stresses (overridden from FEElasticDomain)
	void UpdateStresses(FEModel& fem);

	// update element stress
	void UpdateElementStress(int iel, double dt, bool sstate);

public:
	// internal work (overridden from FEElasticDomain)
	void InternalForces(FEGlobalVector& R);

	//! internal fluid work
	virtual void InternalFluidWork(vector<double>& R, double dt);

	//! internal fluid work (steady-state analysis)
	virtual void InternalFluidWorkSS(vector<double>& R, double dt);

	//! internal solute work
	virtual void InternalSoluteWork(vector<double>& R, double dt);

	//! internal solute work (steady-state analysis)
	virtual void InternalSoluteWorkSS(vector<double>& R, double dt);

public:
	//! calculates the global stiffness matrix for this domain
	virtual void StiffnessMatrix(FESolver* psolver, bool bsymm, const FETimePoint& tp);

	//! calculates the global stiffness matrix for this domain (steady-state case)
	virtual void StiffnessMatrixSS(FESolver* psolver, bool bsymm, const FETimePoint& tp);

protected:
	//! element internal force vector
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

	//! Calculates the internal fluid forces
	bool ElementInternalFluidWork(FESolidElement& elem, vector<double>& fe, double dt);

	//! Calculates the internal fluid forces for steady-state response
	bool ElementInternalFluidWorkSS(FESolidElement& elem, vector<double>& fe, double dt);
	
	//! Calculates the internal solute forces
	bool ElementInternalSoluteWork(FESolidElement& elem, vector<double>& fe, double dt);

	//! Calculates the internal solute forces for steady-state response
	bool ElementInternalSoluteWorkSS(FESolidElement& elem, vector<double>& fe, double dt);
	
	//! calculates the element solute-poroelastic stiffness matrix
	bool ElementBiphasicSoluteStiffness(FESolidElement& el, matrix& ke, bool bsymm, double dt);

	//! calculates the element solute-poroelastic stiffness matrix
	bool ElementBiphasicSoluteStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm, double dt);
	
	//! calculates the solid element stiffness matrix for steady-state response
	void SolidElementStiffness(FESolidElement& el, matrix& ke);
	
	//! material stiffness component
	void ElementBiphasicSoluteMaterialStiffness(FESolidElement& el, matrix& ke);

	//! geometrical stiffness
	void ElementGeometricalStiffness(FESolidElement &el, matrix &ke);

protected: // overridden from FEElasticDomain, but not implemented in this domain
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) {}
	void InertialForces(FEGlobalVector& R, vector<double>& F) {}
	void StiffnessMatrix(FESolver* psolver) {}
	void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) {}
	void MassMatrix(FESolver* psolver, double scale) {}

protected:
	FEBiphasicSolute*	m_pMat;
};
