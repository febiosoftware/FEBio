#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
class FEHolmesMow : public FEElasticMaterial
{
public:
	FEHolmesMow(FEModel* pfem) : FEElasticMaterial(pfem) {}
		
public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio
	double	m_b;	//!< Exponential stiffening coefficient
	double	lam;	//!< first Lame coefficient
	double	mu;		//!< second Lame coefficient
	double	Ha;		//!< aggregate modulus
		
public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;
		
	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;
		
	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	//! data initialization and checking
	bool Validate() override;
		
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
