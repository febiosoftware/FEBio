// ShenoyMaterial.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "FEShenoyMaterial.h"

BEGIN_FECORE_CLASS(FEShenoyMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
	ADD_PARAMETER(m_k , FE_RANGE_GREATER_OR_EQUAL(0.0), "k");
	ADD_PARAMETER(m_Ef, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ef");
	ADD_PARAMETER(m_lamc, FE_RANGE_GREATER_OR_EQUAL(1.0), "lam_c");
	ADD_PARAMETER(m_lamt, FE_RANGE_NOT_EQUAL(0.0), "lam_t");
	ADD_PARAMETER(m_n   , FE_RANGE_GREATER_OR_EQUAL(1.0), "n");
	ADD_PARAMETER(m_m   , FE_RANGE_GREATER_OR_EQUAL(1.0), "m");
END_FECORE_CLASS();

FEShenoyMaterial::FEShenoyMaterial(FEModel* fem) : FEElasticMaterial(fem)
{
	m_mu = 0.0;
	m_k  = 0.0;
	m_Ef = 0.0;
	m_lamc = 1.0;
	m_lamt = 0.0;
	m_n = 5;
	m_m = 1;
}

double FEShenoyMaterial::fiberStress(double lam)
{
	double lam1 = m_lamc - m_lamt * 0.5;
	double lam2 = m_lamc + m_lamt * 0.5;

	if (lam < lam1) return 0.0;
	else if (lam > lam2)
	{
		return m_Ef*((lam2 - lam1) / (m_n + 1.0) + (pow(1.0 + lam - lam2, m_m + 1.0) - 1.0) / (m_m + 1.0));
	}
	else
	{
		return m_Ef*(pow((lam - lam1)/(lam2 - lam1), m_n)*(lam - lam1)) / (m_n + 1);
	}
}

double FEShenoyMaterial::fiberTangent(double lam)
{
	double lam1 = m_lamc - m_lamt * 0.5;
	double lam2 = m_lamc + m_lamt * 0.5;

	if (lam < lam1) return 0.0;
	else if (lam > lam2)
	{
		return m_Ef*pow(1.0 + lam - lam2, m_m);
	}
	else
	{
		return m_Ef*pow((lam - lam1) / (lam2 - lam1), m_n);
	}
}

mat3ds FEShenoyMaterial::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3ds B = pt.LeftCauchyGreen();
	double J = pt.m_J;
	double Jm23 = pow(J, -2.0 / 3.0);

	mat3ds Bbar = B*Jm23;
	mat3dd I(1.0);

	double lam[3] = {0};
	vec3d n[3];
	B.eigen(lam, n);
	lam[0] = (lam[0] >= 0.0 ? sqrt(lam[0]) : 0.0);
	lam[1] = (lam[1] >= 0.0 ? sqrt(lam[1]) : 0.0);
	lam[2] = (lam[2] >= 0.0 ? sqrt(lam[2]) : 0.0);

	// isotropic component of the stress
	mat3ds s_b = Bbar.dev()*(m_mu / J) + I*(m_k * (J - 1.0));

	// anisotropic (fiber) contribution of the stress
	double df[3] = {0};
	df[0] = fiberStress(lam[0]);
	df[1] = fiberStress(lam[1]);
	df[2] = fiberStress(lam[2]);

	mat3ds s_f = (dyad(n[0])*(df[0] * lam[0]) + dyad(n[1])*(df[1] * lam[1]) + dyad(n[2])*(df[2] * lam[2])) / J;

	// add them together
	mat3ds s = s_b + s_f;

	// and done
	return s;
}

tens4ds FEShenoyMaterial::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3ds B = pt.LeftCauchyGreen();
	double J = pt.m_J;
	double Jm23 = pow(J, -2.0 / 3.0);

	mat3ds Bbar = B*Jm23;
	mat3dd I(1.0);

	double lam[3] = { 0 };
	vec3d n[3];
	B.eigen(lam, n);
	lam[0] = (lam[0] >= 0.0 ? sqrt(lam[0]) : 0.0);
	lam[1] = (lam[1] >= 0.0 ? sqrt(lam[1]) : 0.0);
	lam[2] = (lam[2] >= 0.0 ? sqrt(lam[2]) : 0.0);

	double I1bar = Bbar.tr();

	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);

	tens4ds c_b = (I4*(2.0*I1bar / 3.0) - (dyad1s(Bbar, I))*(2.0 / 3.0) + IxI*(I1bar*2.0 / 9.0))*(m_mu / J) \
				+ (IxI*(2.0*J - 1) - I4*(2.0*(J - 1.0)))*m_k;

	double df[3] = {0}, ddf[3] = {0}, sf[3] = {0};
	df[0] = fiberStress(lam[0]);
	df[1] = fiberStress(lam[1]);
	df[2] = fiberStress(lam[2]);

	ddf[0] = fiberTangent(lam[0]);
	ddf[1] = fiberTangent(lam[1]);
	ddf[2] = fiberTangent(lam[2]);

	sf[0] = lam[0] * df[0] / J;
	sf[1] = lam[1] * df[1] / J;
	sf[2] = lam[2] * df[2] / J;

	tens4ds c_f; c_f.zero();
	for (int a=0; a<3; ++a)
	{
		c_f += dyad1s(dyad(n[a]))*(lam[a]*(lam[a]*ddf[a] - df[a])/J);
	}
	
	int p[][2] = {{0,1},{1,2},{0,2}};

	for (int i=0; i < 3; ++i)
	{
		int a = p[i][0];
		int b = p[i][1];

		double g = 0.0;
		if (fabs(lam[a] - lam[b]) <= 1e-10)
		{
			g = 0.5*lam[a]*lam[a]*ddf[a] / J - sf[a];
		}
		else
		{
			g = (sf[a]*lam[b]*lam[b] - sf[b]*lam[a]*lam[a])/(lam[a]*lam[a] - lam[b]*lam[b]);
		}

		mat3ds Nab = dyads(n[a], n[b]);
		c_f += dyad1s(Nab)*(g);
	}

	tens4ds c = c_b + c_f;

	return c;
}
