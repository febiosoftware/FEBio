// FELinearElastic.h: interface for the FELinearElastic class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FELINEARELASTIC_H__D91CFCDE_A6EB_4AF5_B6ED_89A4725528FE__INCLUDED_)
#define AFX_FELINEARELASTIC_H__D91CFCDE_A6EB_4AF5_B6ED_89A4725528FE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! Linear elatic material for small rotations and small deformations 

class FELinearElastic : public FEElasticMaterial
{
public:
	FELinearElastic() {}

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
	DECLARE_REGISTERED(FELinearElastic);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FELINEARELASTIC_H__D91CFCDE_A6EB_4AF5_B6ED_89A4725528FE__INCLUDED_)
