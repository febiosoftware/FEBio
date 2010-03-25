#include "stdafx.h"
#include "FESpectralPermeabilityTransIso.h"


// register the material with the framework
REGISTER_MATERIAL(FESpectralPermeabilityTransIso, "spectral permeability trans iso");

// define the material parameters
BEGIN_PARAMETER_LIST(FESpectralPermeabilityTransIso, FEPoroElastic)
ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
ADD_PARAMETER(m_permA, FE_PARAM_DOUBLE, "permA");
ADD_PARAMETER(m_permT, FE_PARAM_DOUBLE, "permT");
ADD_PARAMETER(m_MA, FE_PARAM_DOUBLE, "MA");
ADD_PARAMETER(m_MT, FE_PARAM_DOUBLE, "MT");
ADD_PARAMETER(m_alphaA, FE_PARAM_DOUBLE, "alphaA");
ADD_PARAMETER(m_alphaT, FE_PARAM_DOUBLE, "alphaT");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FESpectralPermeabilityTransIso::FESpectralPermeabilityTransIso()
{
	m_permA = m_permT = 1;
	m_phi0 = 0;
	m_MA = m_MT = 0;
	m_alphaA = m_alphaT = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FESpectralPermeabilityTransIso::Init()
{
	if (m_permA < 0) throw MaterialError("permA must be >= 0");
	if (m_permT < 0) throw MaterialError("permT must be >= 0");
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 < phi0 <= 1");
	if (m_MA < 0) throw MaterialError("MA must be >= 0");
	if (m_alphaA < 0) throw MaterialError("alphaA must be >= 0");
	if (m_MT < 0) throw MaterialError("MT must be >= 0");
	if (m_alphaT < 0) throw MaterialError("alphaT must be >= 0");
}

//-----------------------------------------------------------------------------
//! Fluid flux.

vec3d FESpectralPermeabilityTransIso::Flux(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d& Q = et.Q;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
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
	double KA, KT;
	vec3d V;
	mat3ds m;
	KA = m_permA*pow((J-m_phi0)/m_phi0,m_alphaA)*exp(m_MA*(J*J-1.0)/2.0);
	KT = m_permT*pow((J-m_phi0)/m_phi0,m_alphaT)*exp(m_MT*(J*J-1.0)/2.0);
	V.x = Q[0][0]; V.y = Q[1][0]; V.z = Q[2][0];
	m = dyad(F*V)/J;
	kt = (KT/J)*b + (KA-KT)*m;
	w = -(kt*gradp);
	
	return w;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
void FESpectralPermeabilityTransIso::Permeability(double k[3][3], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d& Q = et.Q;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
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
	double KA, KT;
	vec3d V;
	mat3ds m;
	KA = m_permA*pow((J-m_phi0)/m_phi0,m_alphaA)*exp(m_MA*(J*J-1.0)/2.0);
	KT = m_permT*pow((J-m_phi0)/m_phi0,m_alphaT)*exp(m_MT*(J*J-1.0)/2.0);
	V.x = Q[0][0]; V.y = Q[1][0]; V.z = Q[2][0];
	m = dyad(F*V)/J;
	kt = (KT/J)*b + (KA-KT)*m;
	k[0][0] = kt.xx();
	k[1][1] = kt.yy();
	k[2][2] = kt.zz();
	k[0][1] = k[1][0] = kt.xy();
	k[1][2] = k[2][1] = kt.yz();
	k[2][0] = k[0][2] = kt.xz();
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4ds FESpectralPermeabilityTransIso::Tangent_Permeability(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d& Q = et.Q;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// deformation gradient
	mat3d &F = et.F;
	
	// relative volume
	double J = et.J;
	
	// --- strain-dependent permeability ---
	
	tens4ds k4;
	mat3dd I(1);	// Identity
	double KA, KT;
	double KAprime, KTprime;
	vec3d V;
	mat3ds m;
	KA = m_permA*pow((J-m_phi0)/m_phi0,m_alphaA)*exp(m_MA*(J*J-1.0)/2.0);
	KT = m_permT*pow((J-m_phi0)/m_phi0,m_alphaT)*exp(m_MT*(J*J-1.0)/2.0);
	KAprime = (m_alphaA/(J-m_phi0) + m_MA*J)*KA;
	KTprime = (m_alphaT/(J-m_phi0) + m_MT*J)*KT;
	V.x = Q[0][0]; V.y = Q[1][0]; V.z = Q[2][0];
	m = dyad(F*V)/J;
	// Divide by two because dyad1s evaluates twice the symmetric part
	k4 = dyad1s(KTprime*b+J*(KAprime-KTprime)*m,I)/2.0;
	
	return k4;
}
