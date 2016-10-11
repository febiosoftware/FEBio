// FETransIsoMooneyRivlin.cpp: implementation of the FETransIsoMooneyRivlin class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FETransIsoMooneyRivlin.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FETransIsoMooneyRivlin, FEUncoupledMaterial)
	ADD_PARAMETER(c1          , FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(c2          , FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETER(m_fib.m_c3  , FE_PARAM_DOUBLE, "c3");
	ADD_PARAMETER(m_fib.m_c4  , FE_PARAM_DOUBLE, "c4");
	ADD_PARAMETER(m_fib.m_c5  , FE_PARAM_DOUBLE, "c5");
	ADD_PARAMETER(m_fib.m_lam1, FE_PARAM_DOUBLE, "lam_max");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FETransIsoMooneyRivlin
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
FETransIsoMooneyRivlin::FETransIsoMooneyRivlin(FEModel* pfem) : FEUncoupledMaterial(pfem), m_fib(pfem)
{
	AddProperty(&m_ac, "active_contraction", false);
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FETransIsoMooneyRivlin::CreateMaterialPointData()
{
	return m_fib.CreateMaterialPointData();
}

//-----------------------------------------------------------------------------
mat3ds FETransIsoMooneyRivlin::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
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
	double W1 = c1;
	double W2 = c2;
	// ------------------------------------------------

	// calculate T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	// calculate stress s = pI + 2/J * dev(T) 
	mat3ds s = T.dev()*(2.0/J);

	// calculate the passive fiber stress
	mat3ds fs = m_fib.DevStress(mp);

	// calculate the active fiber stress
	if ((FEActiveFiberContraction*)m_ac) fs += m_ac->FiberStress(pt);

	return s + fs;
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4ds FETransIsoMooneyRivlin::DevTangent(FEMaterialPoint& mp)
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
	double W1, W2;
	W1 = c1;
	W2 = c2;
	// ------------------------------------

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2;

	// calculate C:d2WdCdC:C
	double CWWC = 2*I2*W2;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// deviatoric cauchy-stress, trs = trace[s]/3
	mat3ds devs = pt.m_s.dev();

	// d2W/dCdC:C
	mat3ds WCCxC = B*(W2*I1) - B2*W2;

	tens4ds cw = (BxB - B4)*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);
	tens4ds c = dyad1s(devs, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	return c + m_fib.DevTangent(mp); // + m_pafc->FiberTangent(mp);
}

//-----------------------------------------------------------------------------
double FETransIsoMooneyRivlin::DevStrainEnergyDensity(FEMaterialPoint& mp)
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
	double sed = c1*(I1-3) + c2*(I2-3);
    
	// add the fiber sed
	sed += m_fib.DevStrainEnergyDensity(mp);
    
	return sed;
}
