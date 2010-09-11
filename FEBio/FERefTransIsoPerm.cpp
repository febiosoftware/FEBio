#include "stdafx.h"
#include "FERefTransIsoPerm.h"


// register the material with the framework
REGISTER_MATERIAL(FERefTransIsoPerm, "ref trans iso perm");

// define the material parameters
BEGIN_PARAMETER_LIST(FERefTransIsoPerm, FEPoroElastic)
ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
ADD_PARAMETER(m_perm0, FE_PARAM_DOUBLE, "perm0");
ADD_PARAMETER(m_perm1T, FE_PARAM_DOUBLE, "perm1T");
ADD_PARAMETER(m_perm1A, FE_PARAM_DOUBLE, "perm1A");
ADD_PARAMETER(m_perm2T, FE_PARAM_DOUBLE, "perm2T");
ADD_PARAMETER(m_perm2A, FE_PARAM_DOUBLE, "perm2A");
ADD_PARAMETER(m_M0, FE_PARAM_DOUBLE, "M0");
ADD_PARAMETER(m_MT, FE_PARAM_DOUBLE, "MT");
ADD_PARAMETER(m_MA, FE_PARAM_DOUBLE, "MA");
ADD_PARAMETER(m_alpha0, FE_PARAM_DOUBLE, "alpha0");
ADD_PARAMETER(m_alphaT, FE_PARAM_DOUBLE, "alphaT");
ADD_PARAMETER(m_alphaA, FE_PARAM_DOUBLE, "alphaA");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FERefTransIsoPerm::FERefTransIsoPerm()
{
	m_perm0 = 1;
	m_perm1T = m_perm1A = 0;
	m_perm2T = m_perm2A = 0;
	m_phi0 = 0.5;
	m_M0 = m_MT = m_MA = 0;
	m_alpha0 = m_alphaT = m_alphaA = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FERefTransIsoPerm::Init()
{
	if (m_perm0 < 0) throw MaterialError("perm0 must be >= 0");
	if (m_perm1T < 0) throw MaterialError("perm1T must be >= 0");
	if (m_perm1A < 0) throw MaterialError("perm1A must be >= 0");
	if (m_perm2T < 0) throw MaterialError("perm2T must be >= 0");
	if (m_perm2A < 0) throw MaterialError("perm2A must be >= 0");
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 < phi0 <= 1");
	if (m_M0 < 0) throw MaterialError("M0 must be >= 0");
	if (m_MT < 0) throw MaterialError("MT must be >= 0");
	if (m_MA < 0) throw MaterialError("MA must be >= 0");
	if (m_alpha0 < 0) throw MaterialError("alpha0 must be >= 0");
	if (m_alphaT < 0) throw MaterialError("alphaT must be >= 0");
	if (m_alphaA < 0) throw MaterialError("alphaA must be >= 0");
}

//-----------------------------------------------------------------------------
//! Fluid flux.

vec3d FERefTransIsoPerm::Flux(FEMaterialPoint& mp)
{
	vec3d V;			// axial material directions in reference configuration
	mat3ds m;			// axial texture tensor in current configuration

	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// deformation gradient
	mat3d &F = et.F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	
	// Copy the texture direction in the reference configuration to V
	V.x = et.Q[0][0]; V.y = et.Q[1][0]; V.z = et.Q[2][0];
	m = dyad(F*V);	// Evaluate texture tensor in the current configuration
	
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();
	
	// pressure gradient
	vec3d gradp = pt.m_gradp;
	
	// fluid flux w = -k*grad(p)
	vec3d w;
	double f, k1T, k1A, k2T, k2A;
	double k0 = m_perm0*pow((J-m_phi0)/(1-m_phi0),m_alpha0)*exp(m_M0*(J*J-1.0)/2.0);
	// Transverse direction
	f = pow((J-m_phi0)/(1-m_phi0),m_alphaT)*exp(m_MT*(J*J-1.0)/2.0);
	k1T = m_perm1T/(J*J)*f;
	k2T = 0.5*m_perm2T/pow(J,4)*f;
	// Axial direction
	f = pow((J-m_phi0)/(1-m_phi0),m_alphaA)*exp(m_MA*(J*J-1.0)/2.0);
	k1A = m_perm1A/(J*J)*f;
	k2A = 0.5*m_perm2A/pow(J,4)*f;
	// flux
	w = (-k0*I
		 -k1T*b - (k1A-k1T)*m
		 -2*k2T*b*b - (k2A-k2T)*(m*b+b*m))*gradp;
	
	return w;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
void FERefTransIsoPerm::Permeability(double k[3][3], FEMaterialPoint& mp)
{
	vec3d V;			// axial material directions in reference configuration
	mat3ds m;			// axial texture tensor in current configuration
	
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// deformation gradient
	mat3d &F = et.F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	
	// Copy the texture direction in the reference configuration to V
	V.x = et.Q[0][0]; V.y = et.Q[1][0]; V.z = et.Q[2][0];
	m = dyad(F*V);	// Evaluate texture tensor in the current configuration
	
	// --- strain-dependent permeability ---
	
	double f, k1T, k1A, k2T, k2A;
	double k0 = m_perm0*pow((J-m_phi0)/(1-m_phi0),m_alpha0)*exp(m_M0*(J*J-1.0)/2.0);
	// Transverse direction
	f = pow((J-m_phi0)/(1-m_phi0),m_alphaT)*exp(m_MT*(J*J-1.0)/2.0);
	k1T = m_perm1T/(J*J)*f;
	k2T = 0.5*m_perm2T/pow(J,4)*f;
	// Axial direction
	f = pow((J-m_phi0)/(1-m_phi0),m_alphaA)*exp(m_MA*(J*J-1.0)/2.0);
	k1A = m_perm1A/(J*J)*f;
	k2A = 0.5*m_perm2A/pow(J,4)*f;
	// Permeability
	mat3ds kt = k0*I + k1T*b + (k1A-k1T)*m + 2*k2T*b*b + (k2A-k2T)*(m*b+b*m);
	k[0][0] = kt.xx();
	k[1][1] = kt.yy();
	k[2][2] = kt.zz();
	k[0][1] = k[1][0] = kt.xy();
	k[1][2] = k[2][1] = kt.yz();
	k[2][0] = k[0][2] = kt.xz();
	
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4ds FERefTransIsoPerm::Tangent_Permeability(FEMaterialPoint &mp)
{
	vec3d V;			// axial material directions in reference configuration
	mat3ds m;			// axial texture tensor in current configuration
	
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// deformation gradient
	mat3d &F = et.F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	
	// Copy the texture direction in the reference configuration to V
	V.x = et.Q[0][0]; V.y = et.Q[1][0]; V.z = et.Q[2][0];
	m = dyad(F*V);	// Evaluate texture tensor in the current configuration
	
	double f, k0, K0prime;
	mat3ds k0hat, k1hat, k2hat;
	k0 = m_perm0*pow((J-m_phi0)/(1-m_phi0),m_alpha0)*exp(m_M0*(J*J-1.0)/2.0);
	K0prime = (J*J*m_M0+(J*(m_alpha0+1)-m_phi0)/(J-m_phi0))*k0;
	k0hat = mat3dd(K0prime);
	tens4ds K4 = dyad1s(I,k0hat)/2.0-dyad4s(I)*2*k0;
	// Transverse direction
	f = pow((J-m_phi0)/(1-m_phi0),m_alphaT)*exp(m_MT*(J*J-1.0)/2.0);
	double k1T = m_perm1T/(J*J)*f;
	double k2T = 0.5*m_perm2T/pow(J,4)*f;
	mat3ds k1hatT = mat3dd((J*J*m_MT+(J*(m_alphaT-1)+m_phi0)/(J-m_phi0))*k1T);
	mat3ds k2hatT = mat3dd((J*J*m_MT+(J*(m_alphaT-3)+3*m_phi0)/(J-m_phi0))*k2T);
	// Axial direction
	f = pow((J-m_phi0)/(1-m_phi0),m_alphaA)*exp(m_MA*(J*J-1.0)/2.0);
	double k1A = m_perm1A/(J*J)*f;
	double k2A = 0.5*m_perm2A/pow(J,4)*f;
	mat3ds k1hatA = mat3dd((J*J*m_MA+(J*(m_alphaA-1)+m_phi0)/(J-m_phi0))*k1A);
	mat3ds k2hatA = mat3dd((J*J*m_MA+(J*(m_alphaA-3)+3*m_phi0)/(J-m_phi0))*k2A);
	//  Tangent
	K4 += dyad1s(b*b,k2hatT) + dyad4s(b)*4*k2T + dyad4s(m,b)*2*(k2A-k2T)
	+ (dyad1s(b,k1hatT) + dyad1s(m,k1hatA-k1hatT) + dyad1s(m*b+b*m,k2hatA-k2hatT))/2.0;
	
	return K4;
}
