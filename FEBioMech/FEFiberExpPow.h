#pragma once
#include "FEElasticFiberMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Exponential-power law

class FEFiberExpPow : public FEElasticFiberMaterial
{
public:
	FEFiberExpPow(FEModel* pfem);
	
	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp) override;
	
	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;
	
	//! Strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp) override;
    
protected:
	double	m_alpha;	// coefficient of (In-1) in exponential
	double	m_beta;		// power of (In-1) in exponential
	double	m_ksi;		// fiber modulus

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Exponential-power law
//! (Variation that includes a shear term)
//! TODO: I want to delete one of these two formulations.
class FEFiberExponentialPower : public FEElasticFiberMaterial
{
public:
	FEFiberExponentialPower(FEModel* pfem);

	//! Initialization
	bool Validate() override;

	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp) override;

	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;

	//! Strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp) override;

public:
	double	m_alpha;	// coefficient of (In-1) in exponential
	double	m_beta;		// power of (In-1) in exponential
	double	m_ksi;		// measure of fiber modulus
	double  m_mu;       // shear modulus

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

