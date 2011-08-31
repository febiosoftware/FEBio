#include "stdafx.h"
#include "FEPermRefOrtho.h"


// register the material with the framework
REGISTER_MATERIAL(FEPermRefOrtho, "perm-ref-ortho");

// define the material parameters
BEGIN_PARAMETER_LIST(FEPermRefOrtho, FEHydraulicPermeability)
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
FEPermRefOrtho::FEPermRefOrtho()
{
	m_perm0 = 1;
	m_perm1[0] = m_perm1[1] = m_perm1[2] = 0;
	m_perm2[0] = m_perm2[1] = m_perm2[2] = 0;
	m_M0 = m_alpha0 = 0;
	m_M[0] = m_M[1] = m_M[2] = 0;
	m_alpha[0] = m_alpha[1] =m_alpha[2] = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FEPermRefOrtho::Init()
{
	if (m_perm0 < 0) throw MaterialError("perm0 must be >= 0");
	if (m_perm1[0] < 0) throw MaterialError("perm1 components must be >= 0");
	if (m_perm1[1] < 0) throw MaterialError("perm1 components must be >= 0");
	if (m_perm1[2] < 0) throw MaterialError("perm1 components must be >= 0");
	if (m_perm2[0] < 0) throw MaterialError("perm2 components must be >= 0");
	if (m_perm2[1] < 0) throw MaterialError("perm2 components must be >= 0");
	if (m_perm2[2] < 0) throw MaterialError("perm2 components must be >= 0");
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
//! Permeability tensor.
mat3ds FEPermRefOrtho::Permeability(FEMaterialPoint& mp)
{
	int a;
	vec3d V;			// orthonormal material directions in reference configuration
	mat3ds m[3];		// texture tensor in current configuration
	
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// deformation gradient
	mat3d &F = et.F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	// referential solid volume fraction
	double phi0 = J*(1-pt.m_phiw);
	
	for (a=0; a<3; a++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to V
		V.x = et.Q[0][a]; V.y = et.Q[1][a]; V.z = et.Q[2][a];
		m[a] = dyad(F*V);	// Evaluate texture tensor in the current configuration
	}
	
	// --- strain-dependent permeability ---
	
	double f, k1[3], k2[3];
	double k0 = m_perm0*pow((J-phi0)/(1-phi0),m_alpha0)*exp(m_M0*(J*J-1.0)/2.0);
	for (a=0; a<3; a++) {
		f = pow((J-phi0)/(1-phi0),m_alpha[a])*exp(m_M[a]*(J*J-1.0)/2.0);
		k1[a] = m_perm1[a]/(J*J)*f;
		k2[a] = 0.5*m_perm2[a]/pow(J,4)*f;
	}
	mat3ds kt = k0*I
	+k1[0]*m[0]+k1[1]*m[1]+k1[2]*m[2]
	+k2[0]*(m[0]*b+b*m[0])+k2[1]*(m[1]*b+b*m[1])+k2[2]*(m[2]*b+b*m[2]);
	
	return kt;
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4ds FEPermRefOrtho::Tangent_Permeability_Strain(FEMaterialPoint &mp)
{
	int a;
	vec3d V;			// orthonormal material directions in reference configuration
	mat3ds m[3];		// texture tensor in current configuration
	
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// deformation gradient
	mat3d &F = et.F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.J;
	// referential solid volume fraction
	double phi0 = J*(1-pt.m_phiw);
	
	for (a=0; a<3; a++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to V
		V.x = et.Q[0][a]; V.y = et.Q[1][a]; V.z = et.Q[2][a];
		m[a] = dyad(F*V);	// Evaluate texture tensor in the current configuration
	}
	
	double f, k0, k1, k2, K0prime, K1prime, K2prime;
	mat3ds k0hat, k1hat, k2hat;
	k0 = m_perm0*pow((J-phi0)/(1-phi0),m_alpha0)*exp(m_M0*(J*J-1.0)/2.0);
	K0prime = (J*J*m_M0+(J*(m_alpha0+1)-phi0)/(J-phi0))*k0;
	k0hat = mat3dd(K0prime);
	tens4ds K4 = dyad1s(I,k0hat)/2.0-dyad4s(I)*2*k0;
	for (a=0; a<3; a++) {
		f = pow((J-phi0)/(1-phi0),m_alpha[a])*exp(m_M[a]*(J*J-1.0)/2.0);
		k1 = m_perm1[a]/(J*J)*f;
		k2 = 0.5*m_perm2[a]/pow(J,4)*f;
		K1prime = (J*J*m_M[a]+(J*(m_alpha[a]-1)+phi0)/(J-phi0))*k1;
		K2prime = (J*J*m_M[a]+(J*(m_alpha[a]-3)+3*phi0)/(J-phi0))*k2;
		k1hat = mat3dd(K1prime);
		k2hat = mat3dd(K2prime);
		K4 += (dyad1s(m[a],k1hat) + dyad1s(m[a]*b+b*m[a],k2hat))/2.0
		+dyad4s(m[a],b)*2.0;
	}
	
	return K4;
}
