#pragma once
#include "FECore/FEElasticMaterial.h"

class FEIsotropicElastic : public FEElasticMaterial
{
public:
	FEIsotropicElastic() {}

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
	DECLARE_REGISTERED(FEIsotropicElastic);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
