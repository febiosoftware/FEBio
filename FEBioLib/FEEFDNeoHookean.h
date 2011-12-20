#pragma once
#include "FECore/FEMaterial.h"
#include "FENeoHookean.h"
#include "FEEllipsoidalFiberDistribution.h"

class FEEFDNeoHookean :	public FEElasticMaterial
{
public:
	FEEFDNeoHookean() {}

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

public:

	FEEllipsoidalFiberDistribution	m_EFD;
	FENeoHookean					m_NH;

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
