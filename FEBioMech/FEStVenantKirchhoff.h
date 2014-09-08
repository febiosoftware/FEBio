// FEStVenantKirchhoff.h: interface for the FEStVenantKirchhoff class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESTVENANTKIRCHHOFF_H__5E5C4041_7BDB_4EE5_B092_8A2E120696AD__INCLUDED_)
#define AFX_FESTVENANTKIRCHHOFF_H__5E5C4041_7BDB_4EE5_B092_8A2E120696AD__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Linear elatic material for large deformations

//! This material can be used when a body undergoes large rotations
//! but small strains.

class FEStVenantKirchhoff : public FEElasticMaterial
{
public:
	FEStVenantKirchhoff(FEModel* pfem) : FEElasticMaterial(pfem) {}

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt);
    
	//! data initialization
	void Init();

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FESTVENANTKIRCHHOFF_H__5E5C4041_7BDB_4EE5_B092_8A2E120696AD__INCLUDED_)
