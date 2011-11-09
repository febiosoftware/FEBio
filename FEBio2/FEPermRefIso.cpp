#include "stdafx.h"
#include "FEPermRefIso.h"


// register the material with the framework
REGISTER_MATERIAL(FEPermRefIso, "perm-ref-iso");

// define the material parameters
BEGIN_PARAMETER_LIST(FEPermRefIso, FEHydraulicPermeability)
	ADD_PARAMETER(m_perm0, FE_PARAM_DOUBLE, "perm0");
	ADD_PARAMETER(m_perm1, FE_PARAM_DOUBLE, "perm1");
	ADD_PARAMETER(m_perm2, FE_PARAM_DOUBLE, "perm2");
	ADD_PARAMETER(m_M, FE_PARAM_DOUBLE, "M");
	ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPermRefIso::FEPermRefIso()
{
	m_perm0 = 1;
	m_perm1 = 0;
	m_perm2 = 0;
	m_phi0 = m_M = m_alpha = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FEPermRefIso::Init()
{
	if (m_perm0 < 0) throw MaterialError("perm0 must be >= 0");
	if (m_perm1 < 0) throw MaterialError("perm1 must be >= 0");
	if (m_perm2 < 0) throw MaterialError("perm2 must be >= 0");
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 < phi0 <= 1");
	if (m_M < 0) throw MaterialError("M must be >= 0");
	if (m_alpha < 0) throw MaterialError("alpha must be >= 0");
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
mat3ds FEPermRefIso::Permeability(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	// referential solid volume fraction
	double phi0 = pt.m_phi0;
	
	// --- strain-dependent permeability ---
	
	double f = pow((J-m_phi0)/(1-phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	double k0 = m_perm0*f;
	double k1 = m_perm1/(J*J)*f;
	double k2 = 0.5*m_perm2/pow(J,4)*f;
	mat3ds kt = k0*I+k1*b+2*k2*b*b;
	
	return kt;
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4ds FEPermRefIso::Tangent_Permeability_Strain(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	// referential solid volume fraction
	double phi0 = pt.m_phi0;
	
	double f = pow((J-phi0)/(1-phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	double k0 = m_perm0*f;
	double k1 = m_perm1/(J*J)*f;
	double k2 = 0.5*m_perm2/pow(J,4)*f;
	double K0prime = (J*J*m_M+(J*(m_alpha+1)-phi0)/(J-phi0))*k0;
	double K1prime = (J*J*m_M+(J*(m_alpha-1)+phi0)/(J-phi0))*k1;
	double K2prime = (J*J*m_M+(J*(m_alpha-3)+3*phi0)/(J-phi0))*k2;
	mat3ds k0hat = I*K0prime;
	mat3ds k1hat = I*K1prime;
	mat3ds k2hat = I*K2prime;
	
	tens4ds K4 = dyad1s(I,k0hat)/2.0-dyad4s(I)*2*k0
	+ dyad1s(b,k1hat)/2.0
	+ dyad1s(b*b,k2hat)+dyad4s(b)*4*k2;
	
	return K4;
}
