#pragma once
#include "FEElasticFiberMaterialUC.h"

//-----------------------------------------------------------------------------
//! Neo-Hookean law

class FEFiberNHUC : public FEElasticFiberMaterialUC
{
public:
	FEFiberNHUC(FEModel* pfem) : FEElasticFiberMaterialUC(pfem) { m_mu = 0; }

	//! Cauchy stress
	mat3ds DevStress(FEMaterialPoint& mp, const vec3d& a0) override;

	// Spatial tangent
	tens4ds DevTangent(FEMaterialPoint& mp, const vec3d& a0) override;

	//! Strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;

public:
	double	m_mu;       // shear modulus

						// declare the parameter list
	DECLARE_FECORE_CLASS();
};

