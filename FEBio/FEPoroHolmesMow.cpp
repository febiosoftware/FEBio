#include "stdafx.h"
#include "FEPoroHolmesMow.h"


// register the material with the framework
REGISTER_MATERIAL(FEPoroHolmesMow, "poroelastic Holmes-Mow");

// define the material parameters
BEGIN_PARAMETER_LIST(FEPoroHolmesMow, FEPoroElastic)
	ADD_PARAMETER(m_perm, FE_PARAM_DOUBLE, "perm");
	ADD_PARAMETER(m_permv[0], FE_PARAM_DOUBLE, "permx");
	ADD_PARAMETER(m_permv[1], FE_PARAM_DOUBLE, "permy");
	ADD_PARAMETER(m_permv[2], FE_PARAM_DOUBLE, "permz");
	ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
	ADD_PARAMETER(m_M, FE_PARAM_DOUBLE, "M");
	ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPoroHolmesMow::FEPoroHolmesMow()
{
	m_perm = 1;
	m_permv[0] = m_permv[1] = m_permv[2] = 1;
	m_phi0 = m_M = m_alpha = 0;
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
	double k = m_perm*pow((J-m_phi0)/(1.0-m_phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	w.x = -k*m_permv[0]*gradp.x;
	w.y = -k*m_permv[1]*gradp.y;
	w.z = -k*m_permv[2]*gradp.z;
	
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
	
	double perm = m_perm*pow((J-m_phi0)/(1.0-m_phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	k[0][0] = perm*m_permv[0];
	k[1][1] = perm*m_permv[1];
	k[2][2] = perm*m_permv[2];
	k[0][1] = k[0][2] = 0;
	k[1][0] = k[1][2] = 0;
	k[2][0] = k[2][1] = 0;
}
