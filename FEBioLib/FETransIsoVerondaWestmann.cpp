// FETransIsoVerondaWestmann.cpp: implementation of the FETransIsoVerondaWestmann class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FETransIsoVerondaWestmann.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FETransIsoVerondaWestmann, FETransverselyIsotropic)
	ADD_PARAMETER(m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c2, FE_PARAM_DOUBLE, "c2");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FETransIsoVerondaWestmann
//////////////////////////////////////////////////////////////////////

mat3ds FETransIsoVerondaWestmann::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.F;
	double J = pt.J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1 = m_c1*m_c2*exp(m_c2*(I1-3));
	double W2 = -0.5*m_c1*m_c2;
	// ------------------------------------------------

	// calculate T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	// calculate stress s = pI + 2/J * dev(T) 
	mat3ds s = T.dev()*(2.0/J);

	// add the fiber stress
	s += m_fib.Stress(mp);

	return s;
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4ds FETransIsoVerondaWestmann::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	const double third = 1.0 / 3.0;

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Ji = 1.0/J;

	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1, W2, W11;
	W1 = m_c1*m_c2*exp(m_c2*(I1-3));
	W2 = -0.5*m_c1*m_c2;
	W11 = m_c2*W1;
	// ---

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2;

	// calculate C:d2WdCdC:C
	double CWWC = W11*I1*I1+2*I2*W2;

	// deviatoric cauchy-stress, trs = trace[s]/3
	mat3ds devs = pt.s.dev();

	mat3dd I(1);	// Identity

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// d2W/dCdC:C
	mat3ds WCCxC = B*(I1*(W11 + W2)) - B2*W2;

	tens4ds cw = BxB*((W11 + W2)*4.0*Ji) - B4*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);

	tens4ds c = dyad1s(devs, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	return c + m_fib.Tangent(mp);
}
