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
	FE3FieldElasticSolidDomain(FEModel* pfem) : FEElasticSolidDomain(pfem) {}

	//! \todo Do I really use this?
	FE3FieldElasticSolidDomain& operator = (FE3FieldElasticSolidDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! initialize class
	bool Initialize(FEModel& fem);

	//! Reset data
	void Reset();

	//! augmentation
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

public: // overridden from FEElasticDomain

	// update stresses
	void UpdateStresses(FEModel& fem);

	// calculate stiffness matrix
	void StiffnessMatrix(FESolver* psolver);

protected:
	//! Dilatational stiffness component for nearly-incompressible materials
	void ElementDilatationalStiffness(FEModel& fem, int iel, matrix& ke);

	//! material stiffness component
	void ElementMaterialStiffness(FEModel& fem, int iel, matrix& ke);

	//! geometrical stiffness (i.e. initial stress)
	void ElementGeometricalStiffness(int iel, matrix& ke);

	//! update the stress of an element
	void UpdateElementStress(int iel);

protected:
	vector<ELEM_DATA>	m_Data;
};
