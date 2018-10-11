#pragma once
#include "FEElasticFiberMaterial.h"

//-----------------------------------------------------------------------------
//! Neo-Hookean law

class FEFiberNH : public FEElasticFiberMaterial
{
public:
	FEFiberNH(FEModel* pfem);

	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp, const vec3d& a0) override;

	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp, const vec3d& a0) override;

	//! Strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;

public:
	double	m_mu;       // shear modulus

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
