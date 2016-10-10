#pragma once
#include "FEElasticFiberMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Power law - linear

class FEFiberPowLinear : public FEElasticFiberMaterial
{
public:
    FEFiberPowLinear(FEModel* pfem);
    
    //! Cauchy stress
    mat3ds Stress(FEMaterialPoint& mp);
    
    // Spatial tangent
    tens4ds Tangent(FEMaterialPoint& mp);
    
    //! Strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp);
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
    
public:
    double	m_E;		// fiber modulus
    double  m_lam0;     // stretch ratio at end of toe region
    double  m_beta;     // power law exponent in toe region
};

//-----------------------------------------------------------------------------
//! Power law toe region - linear
//! TODO: I want to delete one of these materials
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
	double  m_beta;     // power law exponent in toe region

private:
	double  m_I0;       // m_lam0^2
	double  m_ksi;      // power law coefficient in toe region
	double  m_b;        // coefficient in linear region

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
