//
//  FESFDParam.cpp
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 7/16/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#include "stdafx.h"
#include "FESFDSBM.h"
#include "FEMultiphasic.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "FEBioMech/geodesic.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// we store the cos and sin of the angles here
int FESFDSBM::m_nres = 0;
double FESFDSBM::m_cth[NSTH];
double FESFDSBM::m_sth[NSTH];
double FESFDSBM::m_cph[NSTH];
double FESFDSBM::m_sph[NSTH];
double FESFDSBM::m_w[NSTH];

// define the material parameters
BEGIN_FECORE_CLASS(FESFDSBM, FEElasticMaterial)
	ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
	ADD_PARAMETER(m_beta , FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETER(m_ksi0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi0" );
	ADD_PARAMETER(m_rho0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "rho0" );
	ADD_PARAMETER(m_g    , FE_RANGE_GREATER_OR_EQUAL(0.0), "gamma");
	ADD_PARAMETER(m_sbm, "sbm");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// FESphericalFiberDistribution
//-----------------------------------------------------------------------------

bool FESFDSBM::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

	// get the parent material which must be a multiphasic material
	FEMultiphasic* pMP = dynamic_cast<FEMultiphasic*> (GetAncestor());
    if (pMP == 0) return MaterialError("Parent material must be multiphasic");
    
	// extract the local id of the SBM whose density controls Young's modulus from the global id
	m_lsbm = pMP->FindLocalSBMID(m_sbm);
	if (m_lsbm == -1) return MaterialError("Invalid value for sbm");

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FESFDSBM::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// get the element's local coordinate system
	mat3d Q = pt.m_Q;
	
	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, Wl;
	const double eps = 0;
	mat3ds s;
	s.zero();
	
    // calculate material coefficients
    double alpha = m_alpha;
    double beta = m_beta;
	double rhor = spt.m_sbmr[m_lsbm];
    double ksi = FiberModulus(rhor);
    
	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];
		
		// --- quadrant 1,1,1 ---
		
		// rotate to reference configuration
		n0e = Q*n0a;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0)*exp(alpha*pow(In - 1.0, beta));
			
			// calculate the stress
			s += dyad(nt)*(Wl*wn);
		}
		
		// --- quadrant -1,1,1 ---
		n0q = vec3d(-n0a.x, n0a.y, n0a.z);
		
		// rotate to reference configuration
		n0e = Q*n0q;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0)*exp(alpha*pow(In - 1.0, beta));
			
			// calculate the stress
			s += dyad(nt)*(Wl*wn);
		}
		
		// --- quadrant -1,-1,1 ---
		n0q = vec3d(-n0a.x, -n0a.y, n0a.z);
		
		// rotate to reference configuration
		n0e = Q*n0q;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0)*exp(alpha*pow(In - 1.0, beta));
			
			// calculate the stress
			s += dyad(nt)*(Wl*wn);
		}
		
		// --- quadrant 1,-1,1 ---
		n0q = vec3d(n0a.x, -n0a.y, n0a.z);
		
		// rotate to reference configuration
		n0e = Q*n0q;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0)*exp(alpha*pow(In - 1.0, beta));
			
			// calculate the stress
			s += dyad(nt)*(Wl*wn);
		}
	}
	
	// we multiply by two to add contribution from other half-sphere
	return s*(4.0/J);
}

//-----------------------------------------------------------------------------
tens4ds FESFDSBM::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// get the element's local coordinate system
	mat3d Q = pt.m_Q;
	
	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, Wll;
	const double eps = 0;
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	tens4ds N4;
	tens4ds c;
	c.zero();
	
    // calculate material coefficients
    double alpha = m_alpha;
    double beta = m_beta;
	double rhor = spt.m_sbmr[m_lsbm];
    double ksi = FiberModulus(rhor);
    
	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];
		
		// --- quadrant 1,1,1 ---
		
		// rotate to reference configuration
		n0e = Q*n0a;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy derivative
			double pIn = alpha*pow(In - 1.0,beta);
			Wll = beta*ksi*pow(In - 1.0, beta-2.0)*(beta*pIn+beta-1.0)*exp(pIn);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(Wll*wn);
		}
		
		// --- quadrant -1,1,1 ---
		
		n0q = vec3d(-n0a.x, n0a.y, n0a.z);
		
		// rotate to reference configuration
		n0e = Q*n0q;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy derivative
			double pIn = alpha*pow(In - 1.0,beta);
			Wll = beta*ksi*pow(In - 1.0, beta-2.0)*(beta*pIn+beta-1.0)*exp(pIn);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(Wll*wn);
		}
		
		// --- quadrant -1,-1,1 ---
		
		n0q = vec3d(-n0a.x, -n0a.y, n0a.z);
		
		// rotate to reference configuration
		n0e = Q*n0q;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy derivative
			double pIn = alpha*pow(In - 1.0,beta);
			Wll = beta*ksi*pow(In - 1.0, beta-2.0)*(beta*pIn+beta-1.0)*exp(pIn);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(Wll*wn);
		}
		
		// --- quadrant 1,-1,1 ---
		
		n0q = vec3d(n0a.x, -n0a.y, n0a.z);
		
		// rotate to reference configuration
		n0e = Q*n0q;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy derivative
			double pIn = alpha*pow(In - 1.0,beta);
			Wll = beta*ksi*pow(In - 1.0, beta-2.0)*(beta*pIn+beta-1.0)*exp(pIn);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(Wll*wn);
		}
	}
	
	// multiply by two to integrate over other half of sphere
	return c*(2.0*4.0/J);
}

//-----------------------------------------------------------------------------
double FESFDSBM::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	
	// get the element's local coordinate system
	mat3d Q = pt.m_Q;
	
	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, W;
	const double eps = 0;
	double sed = 0.0;
	
    // calculate material coefficients
    double alpha = m_alpha;
    double beta = m_beta;
	double rhor = spt.m_sbmr[m_lsbm];
    double ksi = FiberModulus(rhor);
    
	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];
		
		// --- quadrant 1,1,1 ---
		
		// rotate to reference configuration
		n0e = Q*n0a;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy density
            if (m_alpha > 0)
                W = ksi/m_alpha*(exp(alpha*pow(In - 1.0, beta))-1);
            else
                W = ksi*pow(In - 1.0, beta);
			
			// add to total sed
			sed += W*wn;
		}
		
		// --- quadrant -1,1,1 ---
		n0q = vec3d(-n0a.x, n0a.y, n0a.z);
		
		// rotate to reference configuration
		n0e = Q*n0q;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy density
            if (m_alpha > 0)
                W = ksi/m_alpha*(exp(alpha*pow(In - 1.0, beta))-1);
            else
                W = ksi*pow(In - 1.0, beta);
			
			// add to total sed
			sed += W*wn;
		}
		
		// --- quadrant -1,-1,1 ---
		n0q = vec3d(-n0a.x, -n0a.y, n0a.z);
		
		// rotate to reference configuration
		n0e = Q*n0q;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy density
            if (m_alpha > 0)
                W = ksi/m_alpha*(exp(alpha*pow(In - 1.0, beta))-1);
            else
                W = ksi*pow(In - 1.0, beta);
			
			// add to total sed
			sed += W*wn;
		}
		
		// --- quadrant 1,-1,1 ---
		n0q = vec3d(n0a.x, -n0a.y, n0a.z);
		
		// rotate to reference configuration
		n0e = Q*n0q;
		
		// get the global spatial fiber direction in current configuration
		nt = F*n0e;
		
		// Calculate In = n0e*C*n0e
		In = nt*nt;
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate strain energy density
            if (m_alpha > 0)
                W = ksi/m_alpha*(exp(alpha*pow(In - 1.0, beta))-1);
            else
                W = ksi*pow(In - 1.0, beta);
			
			// add to total sed
			sed += W*wn;
		}
	}
	
	// we multiply by two to add contribution from other half-sphere
	return sed*2.0;
}
