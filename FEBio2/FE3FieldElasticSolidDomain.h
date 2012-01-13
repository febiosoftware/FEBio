#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! The following domain implements the finite element formulation for a three-field
//! volume element. 
class FE3FieldElasticSolidDomain : public FEElasticSolidDomain
{
protected:
	struct ELEM_DATA
	{
		double	eJ;		// average element jacobian
		double	ep;		// average pressure
		double	Lk;		// Lagrangian multiplier
	};

public:
	//! constructor
	FE3FieldElasticSolidDomain(FEMesh* pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat) { m_ntype = FE_3F_SOLID_DOMAIN; }

	//! TODO: do I really use this?
	FE3FieldElasticSolidDomain& operator = (FE3FieldElasticSolidDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! create a clone of this class
	FEDomain* Clone()
	{
		FE3FieldElasticSolidDomain* pd = new FE3FieldElasticSolidDomain(m_pMesh, m_pMat);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh; pd->m_Node = m_Node;
		pd->m_Data = m_Data;
		return pd;
	}

	//! initialize class
	bool Initialize(FEModel& fem);

	// update stresses
	void UpdateStresses(FEModel& fem);

	//! calculates the solid element stiffness matrix
	void ElementStiffness(FEM& fem, int iel, matrix& ke);

	//! augmentation
	bool Augment();

	//! serialize data to archive
	void Serialize(DumpFile& ar);

protected:
	//! Dilatational stiffness component for nearly-incompressible materials
	void DilatationalStiffness(FEM& fem, int iel, matrix& ke);

	//! material stiffness component
	void MaterialStiffness(FEM& fem, int iel, matrix& ke);

	//! geometrical stiffness (i.e. initial stress)
	void GeometricalStiffness(int iel, matrix& ke);

protected:
	vector<ELEM_DATA>	m_Data;
};
