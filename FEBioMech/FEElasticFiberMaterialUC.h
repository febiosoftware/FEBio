#pragma once
#include "FEUncoupledMaterial.h"
#include <FECore/FEVectorGenerator.h>

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterialUC : public FEMaterial
{
public:
    FEElasticFiberMaterialUC(FEModel* pfem);

	// calculate stress in fiber direction a0
	virtual mat3ds DevStress(FEMaterialPoint& mp, const vec3d& a0) = 0;

	// Spatial tangent
	virtual tens4ds DevTangent(FEMaterialPoint& mp, const vec3d& a0) = 0;

	//! Strain energy density
	virtual double DevStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) = 0;
};

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
