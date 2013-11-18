#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
class FEIsotropicElastic : public FEElasticMaterial
{
public:
	FEIsotropicElastic(FEModel* pfem) : FEElasticMaterial(pfem) {}

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

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
