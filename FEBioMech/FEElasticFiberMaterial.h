//
//  FEElasticFiberMaterial.h
//
//  Created by Gerard Ateshian on 11/16/13.
//

#pragma once

#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterial : public FEElasticMaterial
{
public:
    FEElasticFiberMaterial(FEModel* pfem) : FEElasticMaterial(pfem) {}
    
    void SetFiberDirection(FEMaterialPoint& mp, const vec3d n0);
};

//-----------------------------------------------------------------------------
//! Exponential-power law

class FEFiberExponentialPower : public FEElasticFiberMaterial
{
public:
	FEFiberExponentialPower(FEModel* pfem) : FEElasticFiberMaterial(pfem) {
        m_alpha = 0; m_beta = 2; m_ksi = 0; m_mu = 0; }
	
	//! Initialization
	void Init();
	
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
	FEFiberNH(FEModel* pfem) : FEElasticFiberMaterial(pfem) { m_mu = 0; }
	
	//! Initialization
	void Init();
	
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
    FEFiberPowerLinear(FEModel* pfem) : FEElasticFiberMaterial(pfem) {
        m_E = 0; m_lam0 = 1; m_beta = 3; }
    
    //! Initialization
    void Init();
    
    //! Cauchy stress
    mat3ds Stress(FEMaterialPoint& mp);
    
    // Spatial tangent
    tens4ds Tangent(FEMaterialPoint& mp);
    
    //! Strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp);
    
public:
    double	m_E;		// fiber modulus
    double  m_lam0;     // stretch ratio at end of toe region
    double  m_beta;     // power law exponent in toe region
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
};
