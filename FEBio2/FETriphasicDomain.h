#pragma once
#include "FEElasticSolidDomain.h"

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
	void StiffnessMatrix(FENLSolver* psolver);
	
	//! calculates the residual
	void Residual(FENLSolver* psolver, vector<double>& R);
	
	//! initialize class
	bool Initialize(FEModel& fem);
	
	// update stresses
	void UpdateStresses(FEModel& fem);

	//! return element stiffness matrix
	void ElementStiffness(FEModel& fem, int iel, matrix& ke) {
		FESolidElement& el = Element(iel);
		ElementTriphasicStiffness(fem, el, ke);
	}
	
	//! calculate internal equivalent nodal forces
	void InternalForces(FEM& fem, FESolidElement& el, vector<double>& fe);
	
protected:
	//! Calculates the internal fluid forces
	bool ElementInternalFluidWork(FEM& fem, FESolidElement& elem, vector<double>& fe);
	
	//! Calculates the internal fluid forces for steady-state response
	bool ElementInternalFluidWorkSS(FEM& fem, FESolidElement& elem, vector<double>& fe);
	
	//! Calculates the internal solute forces
	bool ElementInternalSoluteWork(FEM& fem, FESolidElement& elem, vector<double>& fe, const int ion);
	
	//! Calculates the internal solute forces for steady-state response
	bool ElementInternalSoluteWorkSS(FEM& fem, FESolidElement& elem, vector<double>& fe, const int ion);
	
	//! calculates the element triphasic stiffness matrix
	bool ElementTriphasicStiffness(FEModel& fem, FESolidElement& el, matrix& ke);
	
	//! calculates the element triphasic stiffness matrix
	bool ElementTriphasicStiffnessSS(FEM& fem, FESolidElement& el, matrix& ke);
	
	//! calculates the solid element stiffness matrix for steady-state response
	void SolidElementStiffness(FEM& fem, FESolidElement& el, matrix& ke);
	
	//! material stiffness component
	void ElementTriphasicMaterialStiffness(FEM& fem, FESolidElement& el, matrix& ke);
	
};
