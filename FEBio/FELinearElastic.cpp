// FELinearElastic.cpp: implementation of the FELinearElastic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FELinearElastic.h"

// register the material with the framework
REGISTER_MATERIAL(FELinearElastic, "linear elastic");

// define the parameter list
BEGIN_PARAMETER_LIST(FELinearElastic, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_PARAM_DOUBLE, "E");
	ADD_PARAMETER(m_v, FE_PARAM_DOUBLE, "v");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FELinearElastic
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void FELinearElastic::Init()
{
	// intialize base class
	FEElasticMaterial::Init();

	if (m_E <= 0) throw MaterialError("Invalid value for E");
	if (!INRANGE(m_v, -1.0, 0.5)) throw MaterialError("Invalid value for v");
}

//-----------------------------------------------------------------------------
mat3ds FELinearElastic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;
	double detF = pt.J;

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	// small strain voigt vector
	mat3ds e;

	// caculate small strain tensor
	e.xx() = F[0][0] - 1.0;
	e.yy() = F[1][1] - 1.0;
	e.zz() = F[2][2] - 1.0;
	e.xy() = 0.5*(F[0][1] + F[1][0]);
	e.xz() = 0.5*(F[0][2] + F[2][0]);
	e.yz() = 0.5*(F[1][2] + F[2][1]);

	// return stress
	return mat3ds(1,1,1,0,0,0)*(lam*e.tr()) + e*(2.0*mu);
}

//-----------------------------------------------------------------------------
tens4ds FELinearElastic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	double D[6][6] = {0};
	D[0][0] = lam+2.*mu; D[0][1] = lam      ; D[0][2] = lam      ;
	D[1][0] = lam      ; D[1][1] = lam+2.*mu; D[1][2] = lam      ;
	D[2][0] = lam      ; D[2][1] = lam      ; D[2][2] = lam+2.*mu;
	D[3][3] = mu;
	D[4][4] = mu;
	D[5][5] = mu;

	return tens4ds(D);
}
