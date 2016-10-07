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
	mat3ds Stress(FEMaterialPoint& mp);
	
	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp);
	
	//! Strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp);
    
protected:
	double	m_alpha;	// coefficient of (In-1) in exponential
	double	m_beta;		// power of (In-1) in exponential
	double	m_ksi;		// fiber modulus

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
