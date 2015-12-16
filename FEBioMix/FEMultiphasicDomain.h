#pragma once
#include "FECore/FESolidDomain.h"
#include "FEMultiphasic.h"
#include "FEBioMech/FEElasticDomain.h"
#include "FECore/FETypes.h"

//-----------------------------------------------------------------------------
//! Domain class for multiphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since this domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEMultiphasicDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FEMultiphasicDomain(FEModel* pfem);
	
	//! Reset data
	void Reset();

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat);

	//! Unpack solid element data (overridden from FEDomain)
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! initialize elements for this domain
	void InitElements();
	
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver, bool bsymm, const FETimePoint& tp);

	//! calculates the global stiffness matrix for this domain (steady-state case)
	void StiffnessMatrixSS(FESolver* psolver, bool bsymm, const FETimePoint& tp);
	
	//! initialize class
	bool Initialize(FEModel& fem);

	//! activate
	void Activate();
	
	// update stresses (overridden from FEElasticDomain)
	void UpdateStresses(FEModel& fem);

	// update element state data
	void UpdateElementStress(int iel, double dt);

public:

	// internal work (overridden from FEElasticDomain)
	void InternalForces(FEGlobalVector& R);

	// internal fluid work
	void InternalFluidWork(vector<double>& R, double dt);

	// internal fluid work (steady state analysis)
	void InternalFluidWorkSS(vector<double>& R, double dt);

	// solute work
	void InternalSoluteWork(vector<double>& R, double dt);

	// solute work (steady state analysis)
	void InternalSoluteWorkSS(vector<double>& R, double dt);

public:
	//! element internal force vector
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

	//! Calculates the internal fluid forces
	bool ElementInternalFluidWork(FESolidElement& elem, vector<double>& fe, double dt);
	
	//! Calculates the internal fluid forces for steady-state response
	bool ElementInternalFluidWorkSS(FESolidElement& elem, vector<double>& fe, double dt);
	
	//! Calculates the internal solute forces
	bool ElementInternalSoluteWork(FESolidElement& elem, vector<double>& fe, double dt, const int ion);
	
	//! Calculates the internal solute forces for steady-state response
	bool ElementInternalSoluteWorkSS(FESolidElement& elem, vector<double>& fe, double dt, const int ion);
	
	//! calculates the element triphasic stiffness matrix
	bool ElementMultiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm, double dt);
	
	//! calculates the element triphasic stiffness matrix
	bool ElementMultiphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm, double dt);
	
	//! calculates the solid element stiffness matrix for steady-state response
	void SolidElementStiffness(FESolidElement& el, matrix& ke);
	
	//! material stiffness component
	void ElementMultiphasicMaterialStiffness(FESolidElement& el, matrix& ke);

	//! geometrical stiffness
	void ElementGeometricalStiffness(FESolidElement &el, matrix &ke);

protected: // overridden from FEElasticDomain, but not implemented in this domain
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) {}
	void InertialForces(FEGlobalVector& R, vector<double>& F) {}
	void StiffnessMatrix(FESolver* psolver) {}
	void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) {}
	void MassMatrix(FESolver* psolver, double scale) {}

protected:
	FEMultiphasic*	m_pMat;
	int				m_dofP;		//!< pressure dof index
	int				m_dofC;		//!< concentration dof index
};
