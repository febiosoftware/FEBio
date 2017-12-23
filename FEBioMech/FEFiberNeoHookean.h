#pragma once
#include "FEElasticFiberMaterial.h"

//-----------------------------------------------------------------------------
//! Neo-Hookean law

class FEFiberNH : public FEElasticFiberMaterial
{
public:
	FEFiberNH(FEModel* pfem);

	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp) override;

	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;

	//! Strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp) override;

public:
	double	m_mu;       // shear modulus

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
