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
	enum { MAX_INT = 45 };

public:
	FEEllipsoidalFiberDistribution(FEModel* pfem) : FEElasticMaterial(pfem) {}

	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp) override;

	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp) override;
	
	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
public:
	FEParamDouble	m_beta[3];	// power in power-law relation
	FEParamDouble	m_ksi[3];	// coefficient in power-law relation

								// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Material class for the ellipsoidal fiber distribution
//! (this is the old, obsolete implementation)

class FEEllipsoidalFiberDistributionOld : public FEElasticMaterial
{
public:
	FEEllipsoidalFiberDistributionOld(FEModel* pfem);
	
	//! Initialization
	bool Init() override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp) override;

	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp) override;
	
	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;

protected:
	void InitIntegrationRule();
    
public: // parameters
	FEParamDouble	m_beta[3];	// power in power-law relation
	FEParamDouble	m_ksi[3];	// coefficient in power-law relation

private:
	int		m_nres;	// integration rule
	double	m_cth[NSTH];
	double	m_sth[NSTH];
	double	m_cph[NSTH];
	double	m_sph[NSTH];
	double	m_w[NSTH];

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
