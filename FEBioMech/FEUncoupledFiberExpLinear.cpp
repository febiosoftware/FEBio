#include "stdafx.h"
#include "FEUncoupledFiberExpLinear.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEUncoupledFiberExpLinear, FEUncoupledMaterial);
	ADD_PARAMETER(m_c3  , FE_PARAM_DOUBLE, "c3");
	ADD_PARAMETER(m_c4  , FE_PARAM_DOUBLE, "c4");
	ADD_PARAMETER(m_c5  , FE_PARAM_DOUBLE, "c5");
	ADD_PARAMETER(m_lam1, FE_PARAM_DOUBLE, "lambda");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEUncoupledFiberExpLinear::FEUncoupledFiberExpLinear(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	m_c3 = m_c4 = m_c5 = 0;
	m_lam1 = 1;
}

//-----------------------------------------------------------------------------
//! Fiber material stress
mat3ds FEUncoupledFiberExpLinear::DevStress(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Ji = 1.0 / J;
	double Jm13 = pow(J, -1.0 / 3.0);
	double twoJi = 2.0*Ji;

	// get the initial fiber direction
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// invariant I4
	double I4 = lamd*lamd;

	// strain energy derivative
	double W4 = 0;
	if (lamd > 1)
	{
		double lamdi = 1.0 / lamd;
		double Wl;
		if (lamd < m_lam1)
		{
			Wl = lamdi*m_c3*(exp(m_c4*(lamd - 1)) - 1);
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1 - 1)) - 1) - m_c5*m_lam1;
			Wl = lamdi*(m_c5*lamd + c6);
		}
		W4 = 0.5*lamdi*Wl;
	}
	else
	{
		W4 = 0;
	}

	// calculate dyad of a: AxA = (a x a)
	mat3ds AxA = dyad(a);

	// ---
	// calculate FdWf/dCFt = I4*W4*(a x a)
	mat3ds T = AxA*(W4*I4);

	// calculate stress: 
	mat3ds s = T.dev()*twoJi;

	return s;
}

//-----------------------------------------------------------------------------
//! Fiber material tangent
tens4ds FEUncoupledFiberExpLinear::DevTangent(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);
	double Ji = 1.0 / J;

	// get initial local material axis
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];

	// calculate current local material axis
	vec3d a = F*a0;

	double lam = a.unit();

	// deviatoric stretch
	double lamd = lam*Jm13;

	double I4 = lamd*lamd;

	double W4, W44;
	if (lamd >= 1)
	{
		double lamdi = 1.0 / lamd;
		double Wl, Wll;
		if (lamd < m_lam1)
		{
			Wl = lamdi*m_c3*(exp(m_c4*(lamd - 1)) - 1);
			Wll = m_c3*lamdi*(m_c4*exp(m_c4*(lamd - 1)) - lamdi*(exp(m_c4*(lamd - 1)) - 1));
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1 - 1)) - 1) - m_c5*m_lam1;
			Wl = lamdi*(m_c5*lamd + c6);
			Wll = -c6*lamdi*lamdi;
		}
		W4 = 0.5*lamdi*Wl;
		W44 = 0.25*lamdi*lamdi*(Wll - lamdi*Wl);
	}
	else
	{
		W4 = 0;
		W44 = 0;
	}

	// --- calculate tangent ---

	// calculate dWdC:C
	double WC = W4*I4;

	// calculate C:d2WdCdC:C
	double CWWC = W44*I4*I4;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds Id4 = dyad4s(I);

	mat3ds AxA = dyad(a);
	tens4ds AxAxAxA = dyad1s(AxA);

	tens4ds cw = AxAxAxA*(4.0*Ji*W44*I4*I4) - dyad1s(I, AxA)*(4.0 / 3.0*Ji*W44*I4*I4);

	tens4ds c = (Id4 - IxI / 3.0)*(4.0 / 3.0*Ji*WC) + IxI*(4.0 / 9.0*Ji*CWWC) + cw;

	return c;
}

//-----------------------------------------------------------------------------
//! Fiber material strain energy density
double FEUncoupledFiberExpLinear::DevStrainEnergyDensity(FEMaterialPoint &mp)
{
	// TODO: implement this
	return 0.0;
}
