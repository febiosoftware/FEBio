// FEVerondaWestmann.cpp: implementation of the FEVerondaWestmann class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEVerondaWestmann.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEVerondaWestmann, FEUncoupledMaterial)
	ADD_PARAMETER(m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c2, FE_PARAM_DOUBLE, "c2");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEVerondaWestmann::Init()
{
	FEUncoupledMaterial::Init();

	if (m_c1 <= 0) throw MaterialError("c1 must be positive.");
	if (m_c2 <= 0) throw MaterialError("c2 must be positive.");
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric stress
mat3ds FEVerondaWestmann::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient and its determinant
	mat3d &F = pt.m_F;
	double J = pt.m_J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// --- TODO: put strain energy derivatives here ---
	//
	// W = C1*(exp(C2*(I1-3)-1)-0.5*C1*C2*(I2 - 3)
	//
	// Wi = dW/dIi
	double W1 = m_c1*m_c2*exp(m_c2*(I1-3));
	double W2 = -0.5*m_c1*m_c2;
	// ---

	// calculate T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	return T.dev()*(2.0/J);
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4ds FEVerondaWestmann::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
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
	mat3ds devs = pt.m_s.dev();

	mat3ds I(1,1,1,0,0,0);	// Identity

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// d2W/dCdC:C
	mat3ds WCCxC = B*(I1*(W11 + W2)) - B2*W2;

	tens4ds cw = BxB*((W11 + W2)*4.0*Ji) - B4*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);

	tens4ds c = dyad1s(devs, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEVerondaWestmann::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	// get the elastic material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();
    
	// calculate square of B
	mat3ds B2 = B*B;
    
	// Invariants of B (= invariants of C)
	double I1 = B.tr();
    double I2 = (I1*I1 - B2.tr())/2.0;
    
    double sed = m_c1*(exp(m_c2*(I1-3))-1) - m_c1*m_c2*(I2-3)/2;
    
    return sed;
}
