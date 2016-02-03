#include "stdafx.h"
#include "FEGentMaterial.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEGentMaterial, FEElasticMaterial)
	ADD_PARAMETER2(m_G, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "G");
	ADD_PARAMETER2(m_K, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "K");
	ADD_PARAMETER2(m_Jm, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "Jm");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEGentMaterial::FEGentMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_G = 0.0;
	m_K = 0.0;
	m_Jm = 0.0;
}

//-----------------------------------------------------------------------------
mat3ds FEGentMaterial::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	double J = ep.m_J;
	mat3ds b = ep.RightCauchyGreen();
	double I1 = b.tr();

	double mu = m_G;
	double k  = m_K;
	double Jm = m_Jm;

	double W1 = 0.5*mu*Jm / (Jm - I1 + 3.0);

	double h = 0.5*(J*J - 1.0) - log(J);
	double WJ = 2.0*k*(h*h*h)*(J - 1.0/J);

	mat3ds I(1.0);
	mat3ds s = b*(2.0*W1/J) + I*WJ;

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEGentMaterial::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	double J = ep.m_J;
	mat3ds b = ep.RightCauchyGreen();
	double I1 = b.tr();

	double mu = m_G;
	double k  = m_K;
	double Jm = m_Jm;

	double W11 = 0.5*mu*Jm/((Jm-I1+3)*(Jm-I1+3));

	double h = 0.5*(J*J - 1.0) - log(J);
	double WJ = 2.0*k*(h*h*h)*(J - 1.0/J);
	double WJJ = 6*k*h*h*(J-1.0/J)*(J-1.0/J) + 2*k*h*h*h*(1.0 + 1.0/(J*J));

	mat3ds I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(b);
	
	tens4ds c = BxB*(4.0*W11/J) + IxI*(WJ + J*WJJ) - I4*(2*WJ);

	return c;
}
