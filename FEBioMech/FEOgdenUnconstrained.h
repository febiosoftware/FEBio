#pragma once
#include "FEElasticMaterial.h"

class FEOgdenUnconstrained : public FEElasticMaterial
{
public:
	enum { MAX_TERMS = 6 };
public:
	FEOgdenUnconstrained(FEModel* pfem);
	
	//! calculate the stress
	mat3ds Stress(FEMaterialPoint& pt) override;
	
	//! calculate the tangent
	tens4ds Tangent(FEMaterialPoint& pt) override;
	
	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
protected:
	void EigenValues(mat3ds& A, double l[3], vec3d r[3], const double eps = 0);
	double	m_eps;
	
public:
	double	m_c[MAX_TERMS];		//!< coefficients mu
	double	m_m[MAX_TERMS];		//!< powers
	double	m_p;				//!< coefficient mu prime
	
	DECLARE_FECORE_CLASS();
};
