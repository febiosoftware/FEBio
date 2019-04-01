#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for the spherical fiber distribution
//!

class FESphericalFiberDistribution : public FEElasticMaterial
{
public:
	FESphericalFiberDistribution(FEModel* pfem);
	
	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp) override;
	
	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp) override;
	
	// Strain energy density
	virtual double StrainEnergyDensity(FEMaterialPoint& mp) override;
	
	// declare the parameter list
	DECLARE_FECORE_CLASS();
	
public:
	double	m_beta;		// power in power-law relation
	double	m_ksi;		// coefficient in power-law relation
	double	m_alpha;	// coefficient of exponential argument
};
