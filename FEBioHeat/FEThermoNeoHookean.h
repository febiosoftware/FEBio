#pragma once
#include "FEThermalElastic.h"

class FEThermoNeoHookean : public FEThermalElastic
{
public:
	//! constructor
	FEThermoNeoHookean(FEModel* pfem);

	// initialization
	void Init();

public:
	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp);

	//! spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp);

	//! spatial thermal tangent
	mat3ds ThermalTangent(FEMaterialPoint& mp);

public:
	double	m_mu;		//!< shear modulus
	double	m_k;		//!< bulk modulus
	double	m_a0;		//!< linear expansion coefficient
	double	m_gamma;	//!< exponent of jacobian in energy term

	DECLARE_PARAMETER_LIST();
};
