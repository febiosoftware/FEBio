#include "stdafx.h"
#include "FESpectralPermeabilityIso.h"


// register the material with the framework
REGISTER_MATERIAL(FESpectralPermeabilityIso, "spectral permeability iso");

// define the material parameters
BEGIN_PARAMETER_LIST(FESpectralPermeabilityIso, FEPoroElastic)
ADD_PARAMETER(m_perm, FE_PARAM_DOUBLE, "perm");
ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
ADD_PARAMETER(m_M, FE_PARAM_DOUBLE, "M");
ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FESpectralPermeabilityIso::FESpectralPermeabilityIso()
{
	m_perm = 1;
	m_phi0 = m_M = m_alpha = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FESpectralPermeabilityIso::Init()
{
	if (m_perm < 0) throw MaterialError("perm must be >= 0");
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 < phi0 <= 1");
	if (m_M < 0) throw MaterialError("M must be >= 0");
	if (m_alpha < 0) throw MaterialError("alpha must be >= 0");
}

//-----------------------------------------------------------------------------
//! Fluid flux.

vec3d FESpectralPermeabilityIso::Flux(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();
	
	// pressure gradient
	vec3d gradp = pt.m_gradp;
	
	// fluid flux w = -k*grad(p)
	vec3d w;
	double K = m_perm*pow((J-m_phi0)/m_phi0,m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	w = -(K/J)*b*gradp;
	
	return w;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
void FESpectralPermeabilityIso::Permeability(double k[3][3], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	
	// --- strain-dependent permeability ---
	
	double K = m_perm*pow((J-m_phi0)/m_phi0,m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	mat3ds kt = K/J*b;
	k[0][0] = kt.xx();
	k[1][1] = kt.yy();
	k[2][2] = kt.zz();
	k[0][1] = k[1][0] = kt.xy();
	k[1][2] = k[2][1] = kt.yz();
	k[2][0] = k[0][2] = kt.xz();
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4ds FESpectralPermeabilityIso::Tangent_Permeability(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	
	mat3dd I(1);	// Identity
	tens4ds bxI = dyad1s(b,I)/2.0;
	
	double K = m_perm*pow((J-m_phi0)/m_phi0,m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	double Kprime = (m_alpha/(J-m_phi0) + m_M*J)*K;
	
	return bxI*Kprime;
}
