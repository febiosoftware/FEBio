//
//  FEElasticFiberMaterialUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEElasticFiberMaterialUC.h"
#include "FEFiberMaterialPoint.h"

//-----------------------------------------------------------------------------
// FEFiberExponentialPower
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberExponentialPowerUC, FEElasticFiberMaterialUC)
    ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
    ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
    ADD_PARAMETER(m_ksi , FE_PARAM_DOUBLE, "ksi" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEFiberExponentialPowerUC::Init()
{
    FEMaterial::Init();
	if (m_ksi < 0) throw MaterialError("ksi must be positive.");
	if (m_beta < 2) throw MaterialError("beta must be >= 2.");
	if (m_alpha < 0) throw MaterialError("alpha must be >= 0.");
}

//-----------------------------------------------------------------------------
void FEElasticFiberMaterialUC::SetFiberDirection(FEMaterialPoint& mp, const vec3d n0)
{
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
    pf.m_n0 = n0;
}
//-----------------------------------------------------------------------------
mat3ds FEFiberExponentialPowerUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
	
	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1, Wl;
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;
	
	// fiber direction in global coordinate system
	n0 = pf.m_n0;
	
	// Calculate In = n0*C*n0
	In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		
		// calculate strain energy derivative
		Wl = m_ksi*pow(In_1, m_beta-1.0)*exp(m_alpha*pow(In_1, m_beta));
		
		// calculate the fiber stress
		s = N*(2.0*Wl/J);
	}
	else
	{
		s.zero();
	}
	
	return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEFiberExponentialPowerUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
	
	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	vec3d n0, nt;
	double In_1, Wl, Wll;
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
    mat3ds s;
	tens4ds c;
	
	// fiber direction in global coordinate system
	n0 = pf.m_n0;
	
	// Calculate In = n0*C*n0
	In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		tens4ds NxN = dyad1s(N);
		
		// calculate strain energy derivative
		Wl = m_ksi*pow(In_1, m_beta-1.0)*exp(m_alpha*pow(In_1, m_beta));
		
		// calculate the fiber stress
		s = N*(2.0*Wl/J);
        
		// calculate strain energy 2nd derivative
		double tmp = m_alpha*pow(In_1, m_beta);
		Wll = m_ksi*pow(In_1, m_beta-2.0)*((tmp+1)*m_beta-1.0)*exp(tmp);
		
		// calculate the fiber tangent
		c = NxN*(4.0*Wll/J);
	}
	else
	{
		c.zero();
	}
	
	// This is the final value of the elasticity tensor
    mat3dd I(1);
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	c += ((I4+IxI/3.0)*s.tr() - dyad1s(I,s))*(2./3.)
	- (ddots(IxI, c)-IxI*(c.tr()/3.))/3.;
    
	return c;
}

//-----------------------------------------------------------------------------
// FEFiberNH
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberNHUC, FEElasticFiberMaterialUC)
    ADD_PARAMETER(m_mu, FE_PARAM_DOUBLE, "mu");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEFiberNHUC::Init()
{
    FEMaterial::Init();
	if (m_mu < 0) throw MaterialError("mu must be positive.");
}

//-----------------------------------------------------------------------------
mat3ds FEFiberNHUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
	
	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1;
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;
	
	// fiber direction in global coordinate system
	n0 = pf.m_n0;
	
	// Calculate In = n0*C*n0
	In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		
		// calculate the fiber stress
		s = N*(m_mu*In_1/J);
	}
	else
	{
		s.zero();
	}
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberNHUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
	
	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1;
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
    mat3ds s;
	tens4ds c;
	
	// fiber direction in global coordinate system
	n0 = pf.m_n0;
	
	// Calculate In = n0*C*n0
	In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		tens4ds NxN = dyad1s(N);
		
		// calculate the fiber stress
		s = N*(m_mu*In_1/J);
        
		// calculate the fiber tangent
		c = NxN*(2*m_mu/J);
	}
	else
	{
		c.zero();
	}
	
	// This is the final value of the elasticity tensor
    mat3dd I(1);
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	c += ((I4+IxI/3.0)*s.tr() - dyad1s(I,s))*(2./3.)
	- (ddots(IxI, c)-IxI*(c.tr()/3.))/3.;

	return c;
}
