#pragma once
#include "FEMaterial.h"

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
	virtual void Tangent(double D[6][6], FEMaterialPoint& pt);

	//! data initialization and checking
	void Init();

	//! return bulk modulus
	double BulkModulus() { return m_E/(3.0*(1.0 - 2.0*m_v));}

	// declare as registered
	DECLARE_REGISTERED(FENeoHookean);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
