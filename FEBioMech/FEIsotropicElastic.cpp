#include "stdafx.h"
#include "FEIsotropicElastic.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEIsotropicElastic, FEElasticMaterial)
	ADD_PARAMETER2(m_E, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER2(m_v, FE_PARAM_DOUBLE, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
mat3ds FEIsotropicElastic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.m_F;
	double Ji = 1.0 / pt.m_J;

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

//-----------------------------------------------------------------------------
tens4ds FEIsotropicElastic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double Ji = 1.0 / pt.m_J;

	// lame parameters
	double lam = Ji*(m_v*m_E/((1+m_v)*(1-2*m_v)));
	double mu  = Ji*(0.5*m_E/(1+m_v));

	// left cauchy-green matrix (i.e. the 'b' matrix)
	mat3ds b = pt.LeftCauchyGreen();

	return dyad1s(b)*lam + dyad4s(b)*(2.0*mu);
}

//-----------------------------------------------------------------------------
double FEIsotropicElastic::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

    mat3ds E = (pt.RightCauchyGreen() - mat3dd(1))/2;
    
    double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

    double trE = E.tr();
    double Enorm = E.norm();
    
    double sed = lam*trE*trE/2 + mu*Enorm*Enorm;
    
	return sed;
}

//-----------------------------------------------------------------------------
mat3ds FEIsotropicElastic::PK2Stress(FEMaterialPoint& pt, const mat3ds E)
{
    // lame parameters
    double lam = (m_v*m_E/((1+m_v)*(1-2*m_v)));
    double mu  = (0.5*m_E/(1+m_v));
    
    // Identity
    mat3dd I(1);
    
    // calculate stress
    mat3ds S = I*(E.tr()*lam) + E*(2*mu);
    
    return S;
}

//-----------------------------------------------------------------------------
tens4ds FEIsotropicElastic::MaterialTangent(FEMaterialPoint& pt, const mat3ds E)
{
    // lame parameters
    double lam = (m_v*m_E/((1+m_v)*(1-2*m_v)));
    double mu  = (0.5*m_E/(1+m_v));
    
    // Identity
    mat3dd I(1);
    
    tens4ds c = dyad1s(I)*lam + dyad4s(I)*(2*mu);
    
    return c;
}
