#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! The following domain implements the finite element formulation for a three-field
//! volume element. 
class FEBIOMECH_API FE3FieldElasticSolidDomain : public FEElasticSolidDomain
{
protected:
	struct ELEM_DATA
	{
		double	eJ;		// average element jacobian at intermediate time
		double	ep;		// average pressure
		double	Lk;		// Lagrangian multiplier
        double  eJt;    // average element jacobian at current time
        double  eJp;    // average element jacobian at previous time
	};

public:
	//! constructor
	FE3FieldElasticSolidDomain(FEModel* pfem) : FEElasticSolidDomain(pfem) {}

	//! \todo Do I really use this?
	FE3FieldElasticSolidDomain& operator = (FE3FieldElasticSolidDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! initialize class
	bool Init() override;

    //! initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
	//! Reset data
	void Reset() override;

	//! augmentation
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;
    
public: // overridden from FEElasticDomain

	// update stresses
	void Update(const FETimeInfo& tp) override;

	// calculate stiffness matrix
	void StiffnessMatrix(FESolver* psolver) override;

protected:
	//! Dilatational stiffness component for nearly-incompressible materials
	void ElementDilatationalStiffness(FEModel& fem, int iel, matrix& ke);

	//! material stiffness component
	void ElementMaterialStiffness(int iel, matrix& ke);

	//! geometrical stiffness (i.e. initial stress)
	void ElementGeometricalStiffness(int iel, matrix& ke);

	//! update the stress of an element
	void UpdateElementStress(int iel, const FETimeInfo& tp) override;

protected:
	vector<ELEM_DATA>	m_Data;
};
