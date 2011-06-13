#pragma once
#include "FEBioLib/FEUncoupledMaterial.h"

class FEOgdenMaterial :	public FEUncoupledMaterial
{
public:
	enum { MAX_TERMS = 6 };
public:
	FEOgdenMaterial();

	//! data initialization and checking
	void Init();

	//! calculate the deviatoric stress
	mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate the deviatoric tangent
	virtual tens4ds DevTangent(FEMaterialPoint& pt);

protected:
	void EigenValues(mat3ds& A, double l[3], vec3d r[3], const double eps = 0);
	double	m_eps;

public:
	double	m_c[MAX_TERMS];	//!< coefficients
	double	m_m[MAX_TERMS];	//!< powers

	DECLARE_REGISTERED(FEOgdenMaterial);

	DECLARE_PARAMETER_LIST();
};
