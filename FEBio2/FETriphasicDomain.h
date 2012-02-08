#pragma once
#include "FEBioLib/FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
class FEM;

//-----------------------------------------------------------------------------
//! Domain class for triphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since this domain
//! also needs to calculate elastic stiffness contributions.
//!
class FETriphasicDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FETriphasicDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_TRIPHASIC_DOMAIN; }
	
	FEDomain* Clone()
	{
		FETriphasicDomain* pd = new FETriphasicDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}
	
	//! initialize elements for this domain
	void InitElements();
	
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FENLSolver* psolver, bool bsymm, double dt);

	//! calculates the global stiffness matrix for this domain (steady-state case)
	void StiffnessMatrixSS(FENLSolver* psolver, bool bsymm, double dt);
	
	//! calculates the residual
//	void Residual(FENLSolver* psolver, vector<double>& R);
	
	//! initialize class
	bool Initialize(FEModel& fem);
	
	// update stresses
	void UpdateStresses(FEModel& fem);

	//! return element stiffness matrix
/*	void ElementStiffness(int iel, matrix& ke) {
		FESolidElement& el = Element(iel);
		ElementTriphasicStiffness(el, ke);
	}
*/
public:
	// internal fluid work
	void InternalFluidWork(FENLSolver* psolver, vector<double>& R, double dt);

	// internal fluid work (steady state analysis)
	void InternalFluidWorkSS(FENLSolver* psolver, vector<double>& R, double dt);

	// solute work
	void InternalSoluteWork(FENLSolver* psolver, vector<double>& R, double dt, const int ion);

	// solute work (steady state analysis)
	void InternalSoluteWorkSS(FENLSolver* psolver, vector<double>& R, double dt, const int ion);

public:

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
