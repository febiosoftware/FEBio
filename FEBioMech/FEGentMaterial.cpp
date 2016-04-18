#include "stdafx.h"
#include "FEGentMaterial.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEGentMaterial, FEUncoupledMaterial)
	ADD_PARAMETER2(m_G, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "G");
	ADD_PARAMETER2(m_Jm, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "Jm");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEGentMaterial::FEGentMaterial(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	m_G  = 0.0;
	m_Jm = 0.0;
}

//-----------------------------------------------------------------------------
mat3ds FEGentMaterial::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	double J  = ep.m_J;
	mat3ds b  = ep.DevLeftCauchyGreen();
	double I1 = b.tr();

	double mu = m_G;
	double Jm = m_Jm;

	double W1 = 0.5*mu*Jm / (Jm - I1 + 3.0);

	mat3ds T = b*W1;

	return T.dev()*(2.0/J);
}

//-----------------------------------------------------------------------------
tens4ds FEGentMaterial::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;
	double Ji = 1.0/J;

	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.DevLeftCauchyGreen();

	// Invariants of B (= invariants of C)
	double I1 = B.tr();

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1 = 0.5*m_G*m_Jm / (m_Jm - I1 + 3.0);
	// ---

	// calculate dWdC:C
	double WC = W1*I1;

	// deviatoric cauchy-stress
	mat3ds T = B*W1;
	T = T.dev()*(2.0/J);

	// Identity tensor
	mat3ds I(1,1,1,0,0,0);

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	tens4ds c = dyad1s(T, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC);

	return c;
}

//=============================================================================
BEGIN_PARAMETER_LIST(FECompressibleGentMaterial, FEElasticMaterial)
	ADD_PARAMETER2(m_G, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "G");
	ADD_PARAMETER2(m_K, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "K");
	ADD_PARAMETER2(m_Jm, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "Jm");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FECompressibleGentMaterial::FECompressibleGentMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_G = 0.0;
	m_K = 0.0;
	m_Jm = 0.0;
}

//-----------------------------------------------------------------------------
mat3ds FECompressibleGentMaterial::Stress(FEMaterialPoint& mp)
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

	mat3dd I(1.0);
	mat3ds s = b*(2.0*W1/J) + I*WJ;

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FECompressibleGentMaterial::Tangent(FEMaterialPoint& mp)
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

	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(b);
	
	tens4ds c = BxB*(4.0*W11/J) + IxI*(WJ + J*WJJ) - I4*(2*WJ);

	return c;
}
