#include "stdafx.h"
#include "FESpectralPermeabilityOrtho.h"


// register the material with the framework
REGISTER_MATERIAL(FESpectralPermeabilityOrtho, "spectral permeability ortho");

// define the material parameters
BEGIN_PARAMETER_LIST(FESpectralPermeabilityOrtho, FEPoroElastic)
ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
ADD_PARAMETERV(m_perm, FE_PARAM_DOUBLEV, 3, "perm");
ADD_PARAMETERV(m_M, FE_PARAM_DOUBLEV, 3, "M");
ADD_PARAMETERV(m_alpha, FE_PARAM_DOUBLEV, 3, "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FESpectralPermeabilityOrtho::FESpectralPermeabilityOrtho()
{
	m_perm[0] = m_perm[1] = m_perm[2] = 1;
	m_phi0 = 0;
	m_M[0] = m_M[1] = m_M[2] = 0;
	m_alpha[0] = m_alpha[1] = m_alpha[1] = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FESpectralPermeabilityOrtho::Init()
{
	if (m_perm < 0) throw MaterialError("perm must be >= 0");
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 < phi0 <= 1");
	if (m_M[0] < 0) throw MaterialError("M1 must be >= 0");
	if (m_M[1] < 0) throw MaterialError("M2 must be >= 0");
	if (m_M[2] < 0) throw MaterialError("M3 must be >= 0");
	if (m_alpha[0] < 0) throw MaterialError("alpha1 must be >= 0");
	if (m_alpha[1] < 0) throw MaterialError("alpha2 must be >= 0");
	if (m_alpha[2] < 0) throw MaterialError("alpha3 must be >= 0");
}

//-----------------------------------------------------------------------------
//! Fluid flux.

vec3d FESpectralPermeabilityOrtho::Flux(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d& Q = et.Q;
	
	// deformation gradient
	mat3d &F = et.F;
	
	// relative volume
	double J = et.J;
	
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();
	
	// pressure gradient
	vec3d gradp = pt.m_gradp;
	
	// fluid flux w = -k*grad(p)
	vec3d w;
	mat3ds kt;
	double K;
	vec3d V;
	mat3ds m;
	kt.zero();
	for (int i=0; i<3; i++) {
		V.x = Q[0][i]; V.y = Q[1][i]; V.z = Q[2][i];
		m = dyad(F*V)/J;
		K = m_perm[i]*pow((J-m_phi0)/m_phi0,m_alpha[i])*exp(m_M[i]*(J*J-1.0)/2.0);
		kt += K*m;
	}
	w = -(kt*gradp);
	
	return w;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
void FESpectralPermeabilityOrtho::Permeability(double k[3][3], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d& Q = et.Q;
	
	// deformation gradient
	mat3d &F = et.F;
	
	// relative volume
	double J = et.J;
	
	// --- strain-dependent permeability ---
	
	mat3ds kt;
	mat3ds m;
	double K;
	vec3d V;
	kt.zero();
	for (int i=0; i<3; i++) {
		V.x = Q[0][i]; V.y = Q[1][i]; V.z = Q[2][i];
		m = dyad(F*V)/J;
		K = m_perm[i]*pow((J-m_phi0)/m_phi0,m_alpha[i])*exp(m_M[i]*(J*J-1.0)/2.0);
		kt += K*m;
	}
	k[0][0] = kt.xx();
	k[1][1] = kt.yy();
	k[2][2] = kt.zz();
	k[0][1] = k[1][0] = kt.xy();
	k[1][2] = k[2][1] = kt.yz();
	k[2][0] = k[0][2] = kt.xz();
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4ds FESpectralPermeabilityOrtho::Tangent_Permeability(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d& Q = et.Q;
	
	// deformation gradient
	mat3d &F = et.F;
	
	// relative volume
	double J = et.J;
	
	// --- strain-dependent permeability ---
	
	tens4ds k4;
	mat3dd I(1);	// Identity
	double K;
	double Kprime;
	vec3d V;
	mat3ds m;
	k4.zero();
	for (int i=0; i<3; i++) {
		V.x = Q[0][i]; V.y = Q[1][i]; V.z = Q[2][i];
		m = dyad(F*V)/J;
		K = m_perm[i]*pow((J-m_phi0)/m_phi0,m_alpha[i])*exp(m_M[i]*(J*J-1.0)/2.0);
		Kprime = (m_alpha[i]/(J-m_phi0) + m_M[i]*J)*K;
		k4 += J*Kprime*dyad1s(m,I);
	}
	// Divide by two because dyad1s evaluates twice the symmetric part
	k4 = k4/2.0;
	
	return k4;
}
