// FELinearElastic.cpp: implementation of the FELinearElastic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FELinearElastic.h"

// define the parameter list
BEGIN_PARAMETER_LIST(FELinearElastic, FEElasticMaterial)
	ADD_PARAMETER2(m_E, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER2(m_v, FE_PARAM_DOUBLE, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FELinearElastic
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
bool FELinearElastic::Init()
{
	// intialize base class
	if (FEElasticMaterial::Init() == false) return false;

    m_lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	m_mu  = 0.5*m_E/(1+m_v);

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FELinearElastic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.m_F;

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
	return mat3ds(1,1,1,0,0,0)*(m_lam*e.tr()) + e*(2.0*m_mu);
}

//-----------------------------------------------------------------------------
tens4ds FELinearElastic::Tangent(FEMaterialPoint& mp)
{
	double D[6][6] = {0};
	D[0][0] = m_lam+2.*m_mu; D[0][1] = m_lam      ; D[0][2] = m_lam      ;
	D[1][0] = m_lam      ; D[1][1] = m_lam+2.*m_mu; D[1][2] = m_lam      ;
	D[2][0] = m_lam      ; D[2][1] = m_lam      ; D[2][2] = m_lam+2.*m_mu;
	D[3][3] = m_mu;
	D[4][4] = m_mu;
	D[5][5] = m_mu;

	return tens4ds(D);
}

//-----------------------------------------------------------------------------
double FELinearElastic::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// deformation gradient
	mat3d &F = pt.m_F;
    
	// small strain voigt vector
	mat3ds e;
    
	// caculate small strain tensor
	e.xx() = F[0][0] - 1.0;
	e.yy() = F[1][1] - 1.0;
	e.zz() = F[2][2] - 1.0;
	e.xy() = 0.5*(F[0][1] + F[1][0]);
	e.xz() = 0.5*(F[0][2] + F[2][0]);
	e.yz() = 0.5*(F[1][2] + F[2][1]);
    
    double tre = e.tr();
    double enorm = e.norm();
    double sed = m_lam*tre*tre/2 + m_mu*enorm*enorm;
    
	return sed;
}
