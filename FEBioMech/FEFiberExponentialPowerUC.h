#pragma once
#include "FEElasticFiberMaterialUC.h"

//-----------------------------------------------------------------------------
//! Exponential-power law

class FEFiberExponentialPowerUC : public FEElasticFiberMaterialUC
{
public:
	FEFiberExponentialPowerUC(FEModel* pfem);

	//! Validation
	bool Validate() override;

	//! Cauchy stress
	mat3ds DevStress(FEMaterialPoint& mp, const vec3d& a0) override;

	// Spatial tangent
	tens4ds DevTangent(FEMaterialPoint& mp, const vec3d& a0) override;

	//! Strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;

public:
	double	m_alpha;	// coefficient of (In-1) in exponential
	double	m_beta;		// power of (In-1) in exponential
	double	m_ksi;		// fiber modulus
	double  m_mu;       // shear modulus

						// declare the parameter list
	DECLARE_FECORE_CLASS();
};
