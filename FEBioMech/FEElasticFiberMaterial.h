#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterial : public FEElasticMaterial
{
public:
    FEElasticFiberMaterial(FEModel* pfem);

	FEMaterialPoint* CreateMaterialPointData();

	vec3d GetFiberVector(FEMaterialPoint& mp);

protected:
	// NOTE: Some fiber materials define a theta, phi parameter to define the fiber vector.
	//       Although this is deprecated, for backward compatibility this was feature was moved here
	double	m_thd;
	double	m_phd;

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Exponential-power law

class FEFiberExponentialPower : public FEElasticFiberMaterial
{
public:
	FEFiberExponentialPower(FEModel* pfem);
	
	//! Initialization
	bool Validate();
	
	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp);
	
	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp);
	
	//! Strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp);
    
public:
	double	m_alpha;	// coefficient of (In-1) in exponential
	double	m_beta;		// power of (In-1) in exponential
	double	m_ksi;		// measure of fiber modulus
    double  m_mu;       // shear modulus

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Neo-Hookean law

class FEFiberNH : public FEElasticFiberMaterial
{
public:
	FEFiberNH(FEModel* pfem);
	
	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp);
	
	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp);
	
	//! Strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp);
    
public:
	double	m_mu;       // shear modulus

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Power law toe region - linear

class FEFiberPowerLinear : public FEElasticFiberMaterial
{
public:
    FEFiberPowerLinear(FEModel* pfem);
    
    //! Initialization
    bool Validate();
    
    //! Cauchy stress
    mat3ds Stress(FEMaterialPoint& mp);
    
    // Spatial tangent
    tens4ds Tangent(FEMaterialPoint& mp);
    
    //! Strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp);
    
public:
    double	m_E;		// fiber modulus
    double  m_lam0;     // stretch ratio at end of toe region
    double  m_I0;       // m_lam0^2
    double  m_ksi;      // power law coefficient in toe region
    double  m_beta;     // power law exponent in toe region
    double  m_b;        // coefficient in linear region
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
};
