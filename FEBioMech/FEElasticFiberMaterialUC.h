#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterialUC : public FEUncoupledMaterial
{
public:
    FEElasticFiberMaterialUC(FEModel* pfem);

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() override;

	vec3d GetFiberVector(FEMaterialPoint& mp);

protected:
	// NOTE: Some fiber materials define a theta, phi parameter to define the fiber vector.
	//       Although this is deprecated, for backward compatibility this was feature was moved here
	double	m_thd;
	double	m_phd;

	DECLARE_FECORE_CLASS();
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
	mat3ds DevStress(FEMaterialPoint& mp) override;
	
	// Spatial tangent
	tens4ds DevTangent(FEMaterialPoint& mp) override;
	
	//! Strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp) override;
    
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
	mat3ds DevStress(FEMaterialPoint& mp) override;
	
	// Spatial tangent
	tens4ds DevTangent(FEMaterialPoint& mp) override;
	
	//! Strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp) override;
    
public:
	double	m_mu;       // shear modulus
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
