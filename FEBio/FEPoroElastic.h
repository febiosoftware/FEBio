// FEPoroElastic.h: interface for the FEPoroElastic class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEPOROELASTIC_H__C67341B3_B080_4E3B_8B34_D73AEF86BB33__INCLUDED_)
#define AFX_FEPOROELASTIC_H__C67341B3_B080_4E3B_8B34_D73AEF86BB33__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for poroelastic materials.

class FEPoroElastic : public FENestedMaterial
{
public:
	FEPoroElastic();

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FEPoroElasticMaterialPoint(m_pBase->CreateMaterialPointData());
	}

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! calculate fluid flux
	virtual vec3d Flux(FEMaterialPoint& pt) = 0;

	//! permeability
	virtual void Permeability(double k[3][3], FEMaterialPoint& pt) = 0;

	//! tangent of permeability
	virtual tens4ds Tangent_Permeability(FEMaterialPoint& mp) = 0;
};

#endif // !defined(AFX_FEPOROELASTIC_H__C67341B3_B080_4E3B_8B34_D73AEF86BB33__INCLUDED_)
