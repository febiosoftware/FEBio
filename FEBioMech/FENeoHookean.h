#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Neo Hookean material

//! Implementation of a neo-Hookean hyperelastic material.
class FENeoHookean : public FEElasticMaterial
{
public:
	FENeoHookean(FEModel* pfem) : FEElasticMaterial(pfem) {}

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
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
