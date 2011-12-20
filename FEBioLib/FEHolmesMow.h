#pragma once
#include "FECore/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
class FEHolmesMow : public FEElasticMaterial
{
public:
	FEHolmesMow() {}
		
public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio
	double	m_b;	//!< Exponential stiffening coefficient
	double	lam;	//!< first Lame coefficient
	double	mu;		//!< second Lame coefficient
	double	Ha;		//!< aggregate modulus
		
public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);
		
	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);
		
	//! data initialization and checking
	void Init();
		
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
