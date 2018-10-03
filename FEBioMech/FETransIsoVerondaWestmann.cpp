// FETransIsoVerondaWestmann.cpp: implementation of the FETransIsoVerondaWestmann class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FETransIsoVerondaWestmann.h"

// define the material parameters
BEGIN_FECORE_CLASS(FETransIsoVerondaWestmann, FEUncoupledMaterial)
	ADD_PARAMETER(      m_c1  , "c1");
	ADD_PARAMETER(      m_c2  , "c2");
	ADD_PARAMETER(m_fib.m_c3  , "c3");
	ADD_PARAMETER(m_fib.m_c4  , "c4");
	ADD_PARAMETER(m_fib.m_c5  , "c5");
	ADD_PARAMETER(m_fib.m_lam1, "lam_max");

	ADD_PROPERTY(m_ac, "active_contraction", FEProperty::Optional);

END_FECORE_CLASS();

//////////////////////////////////////////////////////////////////////
// FETransIsoVerondaWestmann
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
FETransIsoVerondaWestmann::FETransIsoVerondaWestmann(FEModel* pfem) : FEUncoupledMaterial(pfem), m_fib(pfem)
{
	m_ac = 0;
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FETransIsoVerondaWestmann::CreateMaterialPointData()
{
	return m_fib.CreateMaterialPointData();
}

//-----------------------------------------------------------------------------
mat3ds FETransIsoVerondaWestmann::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1 = m_c1*m_c2*exp(m_c2*(I1-3));
	double W2 = -0.5*m_c1*m_c2;
	// ------------------------------------------------

	// calculate T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	// calculate stress s = pI + 2/J * dev(T) 
	mat3ds s = T.dev()*(2.0/J);

	// add the passive fiber stress
	s += m_fib.DevStress(mp);

	// add the active fiber stress
	if ((FEActiveFiberContraction*)m_ac) s += m_ac->FiberStress(mp);

	return s;
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4ds FETransIsoVerondaWestmann::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
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

	mat3dd I(1);	// Identity

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// d2W/dCdC:C
	mat3ds WCCxC = B*(I1*(W11 + W2)) - B2*W2;

	tens4ds cw = BxB*((W11 + W2)*4.0*Ji) - B4*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);

	tens4ds c = dyad1s(devs, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	return c + m_fib.DevTangent(mp);
}

//-----------------------------------------------------------------------------
double FETransIsoVerondaWestmann::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();
    
	// calculate square of B
	mat3ds B2 = B*B;
    
	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());
    
	// calculate sed
    double sed = m_c1*(exp(m_c2*(I1-3))-1) - m_c1*m_c2*(I2-3)/2;
    
	// add the fiber strain energy density
	sed += m_fib.DevStrainEnergyDensity(mp);
    
	return sed;
}
