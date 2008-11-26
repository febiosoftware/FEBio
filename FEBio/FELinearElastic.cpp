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

void FELinearElastic::Init()
{
	// intialize base class
	FEElasticMaterial::Init();

	if (m_E <= 0) throw MaterialError("Invalid value for E");
	if (!INRANGE(m_v, -1.0, 0.5)) throw MaterialError("Invalid value for v");
}

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

	// trace of e
	double tre;

	// caculate small strain tensor
	e.xx() = F[0][0] - 1.0;
	e.yy() = F[1][1] - 1.0;
	e.zz() = F[2][2] - 1.0;
	e.xy() = 0.5*(F[0][1] + F[1][0]);
	e.xz() = 0.5*(F[0][2] + F[2][0]);
	e.yz() = 0.5*(F[1][2] + F[2][1]);

	// calculate trace of e
	tre = e.xx() + e.yy() + e.zz();

	// calculate stress
	mat3ds s;
	s.xx() = lam*tre + 2.0*mu*e.xx();
	s.yy() = lam*tre + 2.0*mu*e.yy();
	s.zz() = lam*tre + 2.0*mu*e.zz();
	s.xy() = 2.0*mu*e.xy();
	s.yz() = 2.0*mu*e.yz();
	s.xz() = 2.0*mu*e.xz();

	return s;
}

/*
mat3ds FELinearElastic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d F(pt.F);
	double detF = pt.J;

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	// Identity tensor
	diag3d I(1.0);

	// we don't use this, but let's calculate it anyway
	m_K   = lam + 2.0/3.0*mu;

	// strain tensor
	mat3ds e = F.sym() - I;

	// stress tensor
	mat3ds s = I*(lam*e.tr()) + e*(2.0*mu);

	// stress voigth vector
	return mat3ds(s.xx(), s.yy(), s.zz(), s.xy(), s.yz(), s.xz());
}
*/


void FELinearElastic::Tangent(double D[6][6], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	D[0][0] = lam+2.*mu; D[0][1] = lam      ; D[0][2] = lam      ;
	D[1][0] = lam      ; D[1][1] = lam+2.*mu; D[1][2] = lam      ;
	D[2][0] = lam      ; D[2][1] = lam      ; D[2][2] = lam+2.*mu;
	D[3][3] = mu;
	D[4][4] = mu;
	D[5][5] = mu;
}
