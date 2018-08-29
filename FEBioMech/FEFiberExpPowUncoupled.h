#pragma once
#include "FEElasticFiberMaterialUC.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Exponential-power law

class FEFiberExpPowUncoupled : public FEElasticFiberMaterialUC
{
public:
	FEFiberExpPowUncoupled(FEModel* pfem);
	
	//! Cauchy stress
	virtual mat3ds DevStress(FEMaterialPoint& mp) override;
	
	// Spatial tangent
	virtual tens4ds DevTangent(FEMaterialPoint& mp) override;
	
	//! Strain energy density
	virtual double DevStrainEnergyDensity(FEMaterialPoint& mp) override;
    
protected:
	double	m_alpha;	// coefficient of (In-1) in exponential
	double	m_beta;		// power of (In-1) in exponential
	double	m_ksi;		// fiber modulus

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
