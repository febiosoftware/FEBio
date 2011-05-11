#include "stdafx.h"
#include "FEPoroHolmesMow.h"


// register the material with the framework
REGISTER_MATERIAL(FEPoroHolmesMow, "poroelastic Holmes-Mow");

// define the material parameters
BEGIN_PARAMETER_LIST(FEPoroHolmesMow, FEPoroElastic)
	ADD_PARAMETER(m_perm, FE_PARAM_DOUBLE, "perm");
	ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
	ADD_PARAMETER(m_M, FE_PARAM_DOUBLE, "M");
	ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPoroHolmesMow::FEPoroHolmesMow()
{
	m_perm = 1;
	m_phi0 = m_M = m_alpha = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FEPoroHolmesMow::Init()
{
	if (m_perm < 0) throw MaterialError("perm must be >= 0");
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 < phi0 <= 1");
	if (m_M < 0) throw MaterialError("M must be >= 0");
	if (m_alpha < 0) throw MaterialError("alpha must be >= 0");
}

//-----------------------------------------------------------------------------
//! Fluid flux.

vec3d FEPoroHolmesMow::Flux(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// relative volume
	double J = et.J;
	
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();
	
	// pressure gradient
	vec3d gradp = pt.m_gradp;
	
	// fluid flux w = -k*grad(p)
	vec3d w;
	double k0 = m_perm*pow((J-m_phi0)/(1.0-m_phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	w.x = -k0*gradp.x;
	w.y = -k0*gradp.y;
	w.z = -k0*gradp.z;
	
	return w;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
void FEPoroHolmesMow::Permeability(double k[3][3], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// relative volume
	double J = et.J;
	
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();
	
	// --- strain-dependent isotropic permeability ---
	
	double k0 = m_perm*pow((J-m_phi0)/(1.0-m_phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	k[0][0] = k0;
	k[1][1] = k0;
	k[2][2] = k0;
	k[0][1] = k[0][2] = 0;
	k[1][0] = k[1][2] = 0;
	k[2][0] = k[2][1] = 0;
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4ds FEPoroHolmesMow::Tangent_Permeability(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// relative volume
	double J = et.J;

	mat3dd I(1);	// Identity

	double k0 = m_perm*pow((J-m_phi0)/(1.0-m_phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	double K0prime = (J*J*m_M+(J*(m_alpha+1)-m_phi0)/(J-m_phi0))*k0;
	mat3ds k0hat = I*K0prime;

	return dyad1s(I,k0hat)/2.0-dyad4s(I)*2*k0;
}
