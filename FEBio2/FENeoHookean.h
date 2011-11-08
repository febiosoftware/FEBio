#pragma once
#include "FECore/FEElasticMaterial.h"

class FENeoHookean : public FEElasticMaterial
{
public:
	FENeoHookean() {}

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization and checking
	void Init();

	// declare as registered
	DECLARE_REGISTERED(FENeoHookean);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
