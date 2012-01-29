#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
class FEM;

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
	
	FEDomain* Clone()
	{
		FEBiphasicSoluteDomain* pd = new FEBiphasicSoluteDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}

	//! initialize elements for this domain
	void InitElements();
	
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FENLSolver* psolver);
	
	//! calculates the residual
	void Residual(FENLSolver* psolver, vector<double>& R);

	//! initialize class
	bool Initialize(FEModel& fem);
	
	// update stresses
	void UpdateStresses(FEModel& fem);
	
protected:
	//! Calculates the internal fluid forces
	bool ElementInternalFluidWork(FEM& fem, FESolidElement& elem, vector<double>& fe);

	//! Calculates the internal fluid forces for steady-state response
	bool ElementInternalFluidWorkSS(FEM& fem, FESolidElement& elem, vector<double>& fe);
	
	//! Calculates the internal solute forces
	bool ElementInternalSoluteWork(FEM& fem, FESolidElement& elem, vector<double>& fe);

	//! Calculates the internal solute forces for steady-state response
	bool ElementInternalSoluteWorkSS(FEM& fem, FESolidElement& elem, vector<double>& fe);
	
	//! calculates the element solute-poroelastic stiffness matrix
	bool ElementBiphasicSoluteStiffness(FEM& fem, FESolidElement& el, matrix& ke);

	//! calculates the element solute-poroelastic stiffness matrix
	bool ElementBiphasicSoluteStiffnessSS(FEM& fem, FESolidElement& el, matrix& ke);
	
	//! calculates the solid element stiffness matrix for steady-state response
	void SolidElementStiffness(FEM& fem, FESolidElement& el, matrix& ke);
	
	//! material stiffness component
	void ElementBiphasicSoluteMaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);
};
