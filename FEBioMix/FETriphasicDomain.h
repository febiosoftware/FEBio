#pragma once
#include "FEBiphasicSoluteDomain.h"

//-----------------------------------------------------------------------------
//! Domain class for triphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since this domain
//! also needs to calculate elastic stiffness contributions.
//!
class FETriphasicDomain : public FEBiphasicSoluteDomain
{
public:
	//! constructor
	FETriphasicDomain(FEMesh* pm, FEMaterial* pmat);
	
	//! reset domain data
	void Reset();

	//! initialize class
	bool Initialize(FEModel& fem);

	//! initialize elements for this domain
	void InitElements();
	
	//! calculates the residual
//	void Residual(FESolver* psolver, vector<double>& R);
	
	// update stresses
	void UpdateStresses(FEModel& fem);

	//! update element state data
	void UpdateElementStress(int iel);

public: // --- overridden from FEBiphasicSoluteDomain ---

	// internal fluid work
	void InternalFluidWork(vector<double>& R, double dt);

	// internal fluid work (steady state analysis)
	void InternalFluidWorkSS(vector<double>& R, double dt);

	// solute work
	void InternalSoluteWork(vector<double>& R, double dt);

	// solute work (steady state analysis)
	void InternalSoluteWorkSS(vector<double>& R, double dt);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver, bool bsymm, const FETimePoint& tp);

	//! calculates the global stiffness matrix for this domain (steady-state case)
	void StiffnessMatrixSS(FESolver* psolver, bool bsymm, const FETimePoint& tp);

protected:

	//! Calculates the internal fluid forces
	bool ElementInternalFluidWork(FESolidElement& elem, vector<double>& fe, double dt);
	
	//! Calculates the internal fluid forces for steady-state response
	bool ElementInternalFluidWorkSS(FESolidElement& elem, vector<double>& fe, double dt);
	
	//! Calculates the internal solute forces
	bool ElementInternalSoluteWork(FESolidElement& elem, vector<double>& fe, double dt, const int ion);
	
	//! Calculates the internal solute forces for steady-state response
	bool ElementInternalSoluteWorkSS(FESolidElement& elem, vector<double>& fe, double dt, const int ion);
	
	//! calculates the element triphasic stiffness matrix
	bool ElementTriphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm, double dt);
	
	//! calculates the element triphasic stiffness matrix
	bool ElementTriphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm, double dt);
	
	//! calculates the solid element stiffness matrix for steady-state response
	void SolidElementStiffness(FESolidElement& el, matrix& ke);
	
	//! material stiffness component
	void ElementTriphasicMaterialStiffness(FESolidElement& el, matrix& ke);
	
};
