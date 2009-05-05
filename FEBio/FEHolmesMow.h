#pragma once
#include "FEMaterial.h"

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
		
	//! return bulk modulus
	double BulkModulus() { return m_E/(3.0*(1.0 - 2.0*m_v));}
		
	// declare as registered
	DECLARE_REGISTERED(FEHolmesMow);
		
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
