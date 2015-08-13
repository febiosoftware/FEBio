#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
class FEOgdenMaterial :	public FEUncoupledMaterial
{
public:
	enum { MAX_TERMS = 6 };
public:
	FEOgdenMaterial(FEModel* pfem);

	//! calculate the deviatoric stress
	mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate the deviatoric tangent
	tens4ds DevTangent(FEMaterialPoint& pt);

	//! calculate the deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& pt);
    
protected:
	void EigenValues(mat3ds& A, double l[3], vec3d r[3], const double eps = 0);
	double	m_eps;

public:
	double	m_c[MAX_TERMS];	//!< coefficients
	double	m_m[MAX_TERMS];	//!< powers

	DECLARE_PARAMETER_LIST();
};
