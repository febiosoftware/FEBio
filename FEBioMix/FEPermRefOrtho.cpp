#include "FEPermRefOrtho.h"


// define the material parameters
BEGIN_PARAMETER_LIST(FEPermRefOrtho, FEHydraulicPermeability)
	ADD_PARAMETER(m_perm0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "perm0");
	ADD_PARAMETER(m_M0    , FE_RANGE_GREATER_OR_EQUAL(0.0), "M0");
	ADD_PARAMETER(m_alpha0, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha0");
	ADD_PARAMETER(m_perm1 , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "perm1");
	ADD_PARAMETER(m_perm2 , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "perm2");
	ADD_PARAMETER(m_M     , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "M");
	ADD_PARAMETER(m_alpha , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPermRefOrtho::FEPermRefOrtho(FEModel* pfem) : FEHydraulicPermeability(pfem)
{
	m_perm0 = 1;
	m_perm1[0] = m_perm1[1] = m_perm1[2] = 0;
	m_perm2[0] = m_perm2[1] = m_perm2[2] = 0;
	m_M0 = m_alpha0 = 0;
	m_M[0] = m_M[1] = m_M[2] = 0;
	m_alpha[0] = m_alpha[1] =m_alpha[2] = 0;
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
	mat3d &F = et.m_F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.m_J;
	// referential solid volume fraction
	double phi0 = pt.m_phi0;
	
	for (a=0; a<3; a++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to V
		V.x = et.m_Q[0][a]; V.y = et.m_Q[1][a]; V.z = et.m_Q[2][a];
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
	mat3d &F = et.m_F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.m_J;
	// referential solid volume fraction
	double phi0 = pt.m_phi0;
	
	for (a=0; a<3; a++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to V
		V.x = et.m_Q[0][a]; V.y = et.m_Q[1][a]; V.z = et.m_Q[2][a];
		m[a] = dyad(F*V);	// Evaluate texture tensor in the current configuration
	}
	
	double f, k0, k1, k2, K0prime, K1prime, K2prime;
	mat3ds k0hat, k1hat, k2hat;
	k0 = m_perm0*pow((J-phi0)/(1-phi0),m_alpha0)*exp(m_M0*(J*J-1.0)/2.0);
	K0prime = (1+J*(m_alpha0/(J-m_phi0)+m_M0*J))*k0;
	k0hat = mat3dd(K0prime);
	tens4ds K4 = dyad1s(I,k0hat)/2.0-dyad4s(I)*(2*k0);
	for (a=0; a<3; a++) {
		f = pow((J-phi0)/(1-phi0),m_alpha[a])*exp(m_M[a]*(J*J-1.0)/2.0);
		k1 = m_perm1[a]/(J*J)*f;
		k2 = 0.5*m_perm2[a]/pow(J,4)*f;
		K1prime = (J*J*m_M[a]+(J*(m_alpha[a]-1)+phi0)/(J-phi0))*k1;
		K2prime = (J*J*m_M[a]+(J*(m_alpha[a]-3)+3*phi0)/(J-phi0))*k2;
		k1hat = mat3dd(K1prime);
		k2hat = mat3dd(K2prime);
		K4 += (dyad1s(m[a],k1hat) + dyad1s(m[a]*b+b*m[a],k2hat))/2.0
		+dyad4s(m[a],b)*(2.0*k2);
	}
	
	return K4;
}
