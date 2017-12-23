#pragma once
#include "FEThermalElastic.h"

class FEThermoNeoHookean : public FEThermalElastic
{
public:
	//! constructor
	FEThermoNeoHookean(FEModel* pfem);

public:
	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp) override;

	//! spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;

	//! spatial thermal tangent
	mat3ds ThermalTangent(FEMaterialPoint& mp) override;

public:
	double	m_mu;		//!< shear modulus
	double	m_k;		//!< bulk modulus
	double	m_a0;		//!< linear expansion coefficient
	double	m_gamma;	//!< exponent of jacobian in energy term

	DECLARE_PARAMETER_LIST();
};
