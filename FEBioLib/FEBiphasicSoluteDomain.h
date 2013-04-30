#pragma once
#include "FEElasticSolidDomain.h"
#include <FECore/FETypes.h>

//-----------------------------------------------------------------------------
//! Domain class for biphasic-solute 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since this domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBiphasicSoluteDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEBiphasicSoluteDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_BIPHASIC_SOLUTE_DOMAIN; }
	
	//! Create shallow copy
	FEDomain* Clone();

	//! initialize elements for this domain
	void InitElements();
	
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FENLSolver* psolver, bool bsymm, const FETimePoint& tp);

	//! calculates the global stiffness matrix for this domain (steady-state case)
	void StiffnessMatrixSS(FENLSolver* psolver, bool bsymm, const FETimePoint& tp);
	
	//! calculates the residual
//	void Residual(FENLSolver* psolver, vector<double>& R);

	//! initialize class
	bool Initialize(FEModel& fem);
	
	// update stresses
	void UpdateStresses(FEModel& fem);

	// update element stress
	void UpdateElementStress(int iel, double dt, bool sstate);

public:
	//! internal fluid work
	void InternalFluidWork(FENLSolver* psolver, vector<double>& R, double dt);

	//! internal fluid work (steady-state analysis)
	void InternalFluidWorkSS(FENLSolver* psolver, vector<double>& R, double dt);

	//! internal solute work
	void InternalSoluteWork(FENLSolver* psolver, vector<double>& R, double dt);

	//! internal solute work (steady-state analysis)
	void InternalSoluteWorkSS(FENLSolver* psolver, vector<double>& R, double dt);

protected:
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
};
