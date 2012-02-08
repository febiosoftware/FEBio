#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! Domain class for biphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since the biphasic domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBiphasicSolidDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEBiphasicSolidDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_BIPHASIC_DOMAIN; }

	//! clone this domain
	FEDomain* Clone();
	
	//! initialize class
	bool Initialize(FEModel& fem);
	
/*	//! return element stiffness matrix
	void ElementStiffness(FEModel& fem, int iel, matrix& ke) {
		FESolidElement& el = Element(iel);
		ElementBiphasicStiffness(fem, el, ke);
	}
*/

public: // overrides from FEElasticDomain

	// update stresses
	void UpdateStresses(FEModel& fem);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FENLSolver* psolver, bool bsymm, double dt);

	//! calculates the global stiffness matrix (steady-state case)
	void StiffnessMatrixSS(FENLSolver* psolver, bool bsymm, double dt);
	
	//! calculates the residual
//	void Residual(FENLSolver* psolver, vector<double>& R);

public: // TODO: The following functions are to replace Residual

	//! internal fluid work
	void InternalFluidWork(FENLSolver* psolver, vector<double>& R, double dt);

	//! internal fluid work (steady state analysis)
	void InternalFluidWorkSS(FENLSolver* psolver, vector<double>& R, double dt);

public:
	
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
	
	//! material stiffness component
	void ElementBiphasicMaterialStiffness(FESolidElement& el, matrix& ke);
};
