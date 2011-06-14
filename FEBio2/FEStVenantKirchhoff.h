// FEStVenantKirchhoff.h: interface for the FEStVenantKirchhoff class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESTVENANTKIRCHHOFF_H__5E5C4041_7BDB_4EE5_B092_8A2E120696AD__INCLUDED_)
#define AFX_FESTVENANTKIRCHHOFF_H__5E5C4041_7BDB_4EE5_B092_8A2E120696AD__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Linear elatic material for large deformations

//! This material can be used when a body undergoes large rotations
//! but small strains.

class FEStVenantKirchhoff : public FEElasticMaterial
{
public:
	FEStVenantKirchhoff() {}

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization
	void Init();

	//! return bulk modulus
	double BulkModulus() { return m_E/(3.0*(1.0 - 2.0*m_v));}

	// declare as registered
	DECLARE_REGISTERED(FEStVenantKirchhoff);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FESTVENANTKIRCHHOFF_H__5E5C4041_7BDB_4EE5_B092_8A2E120696AD__INCLUDED_)
