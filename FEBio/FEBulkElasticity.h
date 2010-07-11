#pragma once
#include "FEMaterial.h"

class FEBulkElasticity : public FEElasticMaterial
{
public:
	//! When used on its own (not in a solid mixture), this materials
	//! is intrinsically unstable
	FEBulkElasticity() {m_unstable = true;}
		
public:
	double	m_K;	//!< bulk modulus
		
public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);
		
	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);
		
	//! data initialization and checking
	void Init();
		
	//! return bulk modulus
	double BulkModulus() { return m_K;}
		
	// declare as registered
	DECLARE_REGISTERED(FEBulkElasticity);
		
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
