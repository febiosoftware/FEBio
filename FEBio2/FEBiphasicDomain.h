#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! Domain class for biphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since the biphasic domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBiphasicDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEBiphasicDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_BIPHASIC_DOMAIN; }
	
	FEDomain* Clone()
	{
		FEBiphasicDomain* pd = new FEBiphasicDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		return pd;
	}
	
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FENLSolver* psolver);
	
	//! calculates the residual
	void Residual(FENLSolver* psolver, vector<double>& R);

	//! initialize class
	bool Initialize(FEModel& fem);
	
	// update stresses
	void UpdateStresses(FEModel& fem);

	//! return element stiffness matrix
	void ElementStiffness(FEM& fem, int iel, matrix& ke) {
		FESolidElement& el = Element(iel);
		ElementBiphasicStiffness(fem, el, ke);
	}
	
	//! calculate internal equivalent nodal forces
	void InternalForces(FEM& fem, FESolidElement& el, vector<double>& fe);
	
protected:
	//! Calculates the internal fluid forces
	bool InternalFluidWork(FEM& fem, FESolidElement& elem, vector<double>& fe);
	
	//! Calculates the internal fluid forces for steady-state response
	bool InternalFluidWorkSS(FEM& fem, FESolidElement& elem, vector<double>& fe);
	
	//! calculates the element biphasic stiffness matrix
	bool ElementBiphasicStiffness(FEM& fem, FESolidElement& el, matrix& ke);
	
	//! calculates the element biphasic stiffness matrix for steady-state response
	bool ElementBiphasicStiffnessSS(FEM& fem, FESolidElement& el, matrix& ke);
	
	//! calculates the solid element stiffness matrix
	void SolidElementStiffness(FESolidElement& el, matrix& ke);
	
	//! material stiffness component
	void BiphasicMaterialStiffness(FESolidElement& el, matrix& ke);
};
