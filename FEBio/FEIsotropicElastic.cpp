#include "stdafx.h"
#include "FEIsotropicElastic.h"

// register the material with the framework
REGISTER_MATERIAL(FEIsotropicElastic, "isotropic elastic");

// define the material parameters
BEGIN_PARAMETER_LIST(FEIsotropicElastic, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_PARAM_DOUBLE, "E");
	ADD_PARAMETER(m_v, FE_PARAM_DOUBLE, "v");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FEIsotropicElastic
//////////////////////////////////////////////////////////////////////

void FEIsotropicElastic::Init()
{
	// intialize base class
	FEElasticMaterial::Init();

	if (m_E <= 0) throw MaterialError("Invalid value for E");
	if (!INRANGE(m_v, -1.0, 0.5)) throw MaterialError("Invalid value for v");
}

mat3ds FEIsotropicElastic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.F;
	double Ji = 1.0 / pt.J;

	double trE;

	// lame parameters
	double lam = Ji*(m_v*m_E/((1+m_v)*(1-2*m_v)));
	double mu  = Ji*(0.5*m_E/(1+m_v));

	// calculate left Cauchy-Green tensor (ie. b-matrix)
	mat3ds b = pt.LeftCauchyGreen();

	// calculate trace of Green-Lagrance strain tensor
	trE = 0.5*(b.tr()-3);

	// calculate square of b-matrix
	// (we commented out the matrix components we do not need)
	mat3ds b2 = b*b;

	// calculate stress
	mat3ds s = b*(lam*trE - mu) + b2*mu;

	return s;
}

tens4ds FEIsotropicElastic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.F;
	double Ji = 1.0 / pt.J;

	// lame parameters
	double lam = Ji*(m_v*m_E/((1+m_v)*(1-2*m_v)));
	double mu  = Ji*(0.5*m_E/(1+m_v));

	// left cauchy-green matrix (i.e. the 'b' matrix)
	mat3ds b = pt.LeftCauchyGreen();

	return dyad1s(b)*lam + dyad4s(b)*(2.0*mu);
}
