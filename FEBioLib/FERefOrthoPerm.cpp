#include "stdafx.h"
#include "FERefOrthoPerm.h"


// register the material with the framework
REGISTER_MATERIAL(FERefOrthoPerm, "ref ortho perm");

// define the material parameters
BEGIN_PARAMETER_LIST(FERefOrthoPerm, FEPoroElastic)
ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
ADD_PARAMETER(m_perm0, FE_PARAM_DOUBLE, "perm0");
ADD_PARAMETER(m_M0, FE_PARAM_DOUBLE, "M0");
ADD_PARAMETER(m_alpha0, FE_PARAM_DOUBLE, "alpha0");
ADD_PARAMETERV(m_perm1, FE_PARAM_DOUBLEV, 3, "perm1");
ADD_PARAMETERV(m_perm2, FE_PARAM_DOUBLEV, 3, "perm2");
ADD_PARAMETERV(m_M, FE_PARAM_DOUBLEV, 3, "M");
ADD_PARAMETERV(m_alpha, FE_PARAM_DOUBLEV, 3, "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FERefOrthoPerm::FERefOrthoPerm()
{
	m_perm0 = 1;
	m_perm1[0] = m_perm1[1] = m_perm1[2] = 0;
	m_perm2[0] = m_perm2[1] = m_perm2[2] = 0;
	m_phi0 = 0.5;
	m_M0 = m_alpha0 = 0;
	m_M[0] = m_M[1] = m_M[2] = 0;
	m_alpha[0] = m_alpha[1] =m_alpha[2] = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FERefOrthoPerm::Init()
{
	if (m_perm0 < 0) throw MaterialError("perm0 must be >= 0");
	if (m_perm1[0] < 0) throw MaterialError("perm1 components must be >= 0");
	if (m_perm1[1] < 0) throw MaterialError("perm1 components must be >= 0");
	if (m_perm1[2] < 0) throw MaterialError("perm1 components must be >= 0");
	if (m_perm2[0] < 0) throw MaterialError("perm2 components must be >= 0");
	if (m_perm2[1] < 0) throw MaterialError("perm2 components must be >= 0");
	if (m_perm2[2] < 0) throw MaterialError("perm2 components must be >= 0");
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 < phi0 <= 1");
	if (m_M0 < 0) throw MaterialError("M0 must be >= 0");
	if (m_M[0] < 0) throw MaterialError("M components must be >= 0");
	if (m_M[1] < 0) throw MaterialError("M components must be >= 0");
	if (m_M[2] < 0) throw MaterialError("M components must be >= 0");
	if (m_alpha0 < 0) throw MaterialError("alpha0 must be >= 0");
	if (m_alpha[0] < 0) throw MaterialError("alpha components must be >= 0");
	if (m_alpha[1] < 0) throw MaterialError("alpha components must be >= 0");
	if (m_alpha[2] < 0) throw MaterialError("alpha components must be >= 0");
}

//-----------------------------------------------------------------------------
//! Fluid flux.

vec3d FERefOrthoPerm::Flux(FEMaterialPoint& mp)
{
	int a;
	vec3d V;			// orthonormal material directions in reference configuration
	mat3ds m[3];		// texture tensor in current configuration

	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// deformation gradient
	mat3d &F = et.F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	
	for (a=0; a<3; a++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to V
		V.x = et.Q[0][a]; V.y = et.Q[1][a]; V.z = et.Q[2][a];
		m[a] = dyad(F*V);	// Evaluate texture tensor in the current configuration
	}
	
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();
	
	// pressure gradient
	vec3d gradp = pt.m_gradp;
	
	// fluid flux w = -k*grad(p)
	vec3d w;
	double f, k1[3], k2[3];
	double k0 = m_perm0*pow((J-m_phi0)/(1-m_phi0),m_alpha0)*exp(m_M0*(J*J-1.0)/2.0);
	for (a=0; a<3; a++) {
		f = pow((J-m_phi0)/(1-m_phi0),m_alpha[a])*exp(m_M[a]*(J*J-1.0)/2.0);
		k1[a] = m_perm1[a]/(J*J)*f;
		k2[a] = 0.5*m_perm2[a]/pow(J,4)*f;
	}
	w = (-k0*I
		 -k1[0]*m[0]-k1[1]*m[1]-k1[2]*m[2]
		 -k2[0]*(m[0]*b+b*m[0])-k2[1]*(m[1]*b+b*m[1])-k2[2]*(m[2]*b+b*m[2]))*gradp;
	
	return w;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
void FERefOrthoPerm::Permeability(double k[3][3], FEMaterialPoint& mp)
{
	int a;
	vec3d V;			// orthonormal material directions in reference configuration
	mat3ds m[3];		// texture tensor in current configuration
	
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// deformation gradient
	mat3d &F = et.F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	
	for (a=0; a<3; a++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to V
		V.x = et.Q[0][a]; V.y = et.Q[1][a]; V.z = et.Q[2][a];
		m[a] = dyad(F*V);	// Evaluate texture tensor in the current configuration
	}
	
	// --- strain-dependent permeability ---
	
	double f, k1[3], k2[3];
	double k0 = m_perm0*pow((J-m_phi0)/(1-m_phi0),m_alpha0)*exp(m_M0*(J*J-1.0)/2.0);
	for (a=0; a<3; a++) {
		f = pow((J-m_phi0)/(1-m_phi0),m_alpha[a])*exp(m_M[a]*(J*J-1.0)/2.0);
		k1[a] = m_perm1[a]/(J*J)*f;
		k2[a] = 0.5*m_perm2[a]/pow(J,4)*f;
	}
	mat3ds kt = k0*I
	+k1[0]*m[0]+k1[1]*m[1]+k1[2]*m[2]
	+k2[0]*(m[0]*b+b*m[0])+k2[1]*(m[1]*b+b*m[1])+k2[2]*(m[2]*b+b*m[2]);
	k[0][0] = kt.xx();
	k[1][1] = kt.yy();
	k[2][2] = kt.zz();
	k[0][1] = k[1][0] = kt.xy();
	k[1][2] = k[2][1] = kt.yz();
	k[2][0] = k[0][2] = kt.xz();
	
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4ds FERefOrthoPerm::Tangent_Permeability(FEMaterialPoint &mp)
{
	int a;
	vec3d V;			// orthonormal material directions in reference configuration
	mat3ds m[3];		// texture tensor in current configuration
	
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// deformation gradient
	mat3d &F = et.F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	
	for (a=0; a<3; a++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to V
		V.x = et.Q[0][a]; V.y = et.Q[1][a]; V.z = et.Q[2][a];
		m[a] = dyad(F*V);	// Evaluate texture tensor in the current configuration
	}
	
	double f, k0, k1, k2, K0prime, K1prime, K2prime;
	mat3ds k0hat, k1hat, k2hat;
	k0 = m_perm0*pow((J-m_phi0)/(1-m_phi0),m_alpha0)*exp(m_M0*(J*J-1.0)/2.0);
	K0prime = (1+J*(m_alpha0/(J-m_phi0)+m_M0*J))*k0;
	k0hat = mat3dd(K0prime);
	tens4ds K4 = dyad1s(I,k0hat)/2.0-dyad4s(I)*2*k0;
	for (a=0; a<3; a++) {
		f = pow((J-m_phi0)/(1-m_phi0),m_alpha[a])*exp(m_M[a]*(J*J-1.0)/2.0);
		k1 = m_perm1[a]/(J*J)*f;
		k2 = 0.5*m_perm2[a]/pow(J,4)*f;
		K1prime = (J*J*m_M[a]+(J*(m_alpha[a]-1)+m_phi0)/(J-m_phi0))*k1;
		K2prime = (J*J*m_M[a]+(J*(m_alpha[a]-3)+3*m_phi0)/(J-m_phi0))*k2;
		k1hat = mat3dd(K1prime);
		k2hat = mat3dd(K2prime);
		K4 += (dyad1s(m[a],k1hat) + dyad1s(m[a]*b+b*m[a],k2hat))/2.0
			+dyad4s(m[a],b)*(2.0*k2);
	}
	
	return K4;
}
