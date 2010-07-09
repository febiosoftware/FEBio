#pragma once
#include "FEMaterial.h"
#include "FENeoHookean.h"
#include "FEEllipsoidalFiberDistribution.h"

class FERandomFiberNeoHookean :	public FEElasticMaterial
{
public:
	FERandomFiberNeoHookean() {}

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio
	double	m_ksi[3];	//!< ksi
	double	m_beta[3];	//!< beta

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization and checking
	void Init();

	//! return bulk modulus
	double BulkModulus();

public:

	FEEllipsoidalFiberDistribution	m_EFD;
	FENeoHookean					m_NH;

	// declare as registered
	DECLARE_REGISTERED(FERandomFiberNeoHookean);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
