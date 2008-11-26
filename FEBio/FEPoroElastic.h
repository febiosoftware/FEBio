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
//! Poroelastic material

class FEPoroElastic : public FEMaterial
{
public:
	FEPoroElastic();

	double	m_perm;			//!< permeability
	double	m_permv[3];		//!< permeability for diagonal tensor
	int		m_nSolidMat;	//!< ID of solid material

	FEElasticMaterial*	m_psmat;	// pointer to the solid material class

	// returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FEPoroElasticMaterialPoint(m_psmat->CreateMaterialPointData());
	}

	//! return solid component's bulk modulus
	double BulkModulus() { return m_psmat->BulkModulus(); }

	//! return solid component's density
	double Density () { return m_psmat->Density(); }

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual void Tangent(double D[6][6], FEMaterialPoint& pt);

	//! calculate fluid flux
	virtual vec3d Flux(FEMaterialPoint& pt);

	//! permeability
	virtual void Permeability(double k[3][3], FEMaterialPoint& pt);

	// declare as registered
	DECLARE_REGISTERED(FEPoroElastic);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};


#endif // !defined(AFX_FEPOROELASTIC_H__C67341B3_B080_4E3B_8B34_D73AEF86BB33__INCLUDED_)
