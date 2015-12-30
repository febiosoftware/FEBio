#pragma once
#include "FEElasticMaterial.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

//-----------------------------------------------------------------------------
//! Material class for the ellipsoidal fiber distribution
//!

class FEEllipsoidalFiberDistribution : public FEElasticMaterial
{
public:
	FEEllipsoidalFiberDistribution(FEModel* pfem) : FEElasticMaterial(pfem) {}
	
	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp);

	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp);
	
	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt);
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
	
public:
	double	m_beta[3];	// power in power-law relation
	double	m_ksi[3];	// coefficient in power-law relation
};


//-----------------------------------------------------------------------------
//! Material class for the ellipsoidal fiber distribution
//! (this is the old, obsolete implementation)

class FEEllipsoidalFiberDistributionOld : public FEElasticMaterial
{
public:
	FEEllipsoidalFiberDistributionOld(FEModel* pfem) : FEElasticMaterial(pfem) { m_nres = 0; }
	
	//! Initialization
	bool Init();

	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp);

	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp);
	
	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt);
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
	
public:
	double	m_beta[3];	// power in power-law relation
	double	m_ksi[3];	// coefficient in power-law relation

	int		m_nres;	// integration rule
	double	m_cth[NSTH];
	double	m_sth[NSTH];
	double	m_cph[NSTH];
	double	m_sph[NSTH];
	double	m_w[NSTH];
    bool    m_bfirst;
};
