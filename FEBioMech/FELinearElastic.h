// FELinearElastic.h: interface for the FELinearElastic class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FELINEARELASTIC_H__D91CFCDE_A6EB_4AF5_B6ED_89A4725528FE__INCLUDED_)
#define AFX_FELINEARELASTIC_H__D91CFCDE_A6EB_4AF5_B6ED_89A4725528FE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Linear elatic material for small rotations and small deformations 

class FELinearElastic : public FEElasticMaterial
{
public:
	FELinearElastic(FEModel* pfem) : FEElasticMaterial(pfem) {}

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio
    double  m_lam;  //!< first Lamé constant
    double  m_mu;   //!< second Lamé constant

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt);
    
	//! data initialization
	bool Init();

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FELINEARELASTIC_H__D91CFCDE_A6EB_4AF5_B6ED_89A4725528FE__INCLUDED_)
