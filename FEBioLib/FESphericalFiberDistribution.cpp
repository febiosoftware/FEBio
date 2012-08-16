/*
 *  FESphericalFiberDistribution.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 12/26/11.
 *  Copyright 2011 Columbia University. All rights reserved.
 *
 */

#include "stdafx.h"
#include "FESphericalFiberDistribution.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// we store the cos and sin of the angles here
int FESphericalFiberDistribution::m_nres = 0;
double FESphericalFiberDistribution::m_cth[NSTH];
double FESphericalFiberDistribution::m_sth[NSTH];
double FESphericalFiberDistribution::m_cph[NSTH];
double FESphericalFiberDistribution::m_sph[NSTH];
double FESphericalFiberDistribution::m_w[NSTH];

// define the material parameters
BEGIN_PARAMETER_LIST(FESphericalFiberDistribution, FEElasticMaterial)
ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
ADD_PARAMETER(m_ksi , FE_PARAM_DOUBLE, "ksi" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FESphericalFiberDistribution
//-----------------------------------------------------------------------------

void FESphericalFiberDistribution::Init()
{
	FEElasticMaterial::Init();

	if (m_unstable) throw MaterialError("This fibrous material is unstable (collapses on itself) when used alone.  Combine it in a solid mixture with a material that can serve as a ground matrix.");
	if (m_ksi < 0) throw MaterialError("ksi must be positive.");
	if (m_beta < 2) throw MaterialError("beta must be greater than 2.");
	
	static bool bfirst = true;
	
	if (bfirst)
	{
		// select the integration rule
		const int nint = (m_nres == 0? NSTL  : NSTH  );
		const double* phi = (m_nres == 0? PHIL  : PHIH  );
		const double* the = (m_nres == 0? THETAL: THETAH);
		const double* w   = (m_nres == 0? AREAL : AREAH );
		
		for (int n=0; n<nint; ++n)
		{
			m_cth[n] = cos(the[n]);
			m_sth[n] = sin(the[n]);
			m_cph[n] = cos(phi[n]);
			m_sph[n] = sin(phi[n]);
			m_w[n] = w[n];
		}
		
		bfirst = false;
	}
}

//-----------------------------------------------------------------------------
mat3ds FESphericalFiberDistribution::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	
	// get the element's local coordinate system
	mat3d QT = (pt.Q).transpose();
	
	// loop over all integration points
	vec3d n0e, n0a, nt;
	double In, Wl;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds s;
	s.zero();
	
	const int nint = (m_nres == 0? NSTL  : NSTH  );
	
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in reference configuration
		n0e.x = m_cth[n]*m_sph[n];
		n0e.y = m_sth[n]*m_sph[n];
		n0e.z = m_cph[n];
		
		// Calculate In = n0e*C*n0e
		In = n0e*(C*n0e);
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// get the global spatial fiber direction in current configuration
			nt = (F*n0e)/sqrt(In);
			
			// calculate the outer product of nt
			mat3ds N = dyad(nt);
			
			// get the local material fiber direction in reference configuration
			n0a = QT*n0e;
			
			// calculate strain energy derivative
			Wl = m_beta*m_ksi*pow(In - 1.0, m_beta-1.0);
			
			// calculate the stress
			s += N*(2.0/J*In*Wl*m_w[n]);
		}
		
	}
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FESphericalFiberDistribution::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	
	// get the element's local coordinate system
	mat3d QT = (pt.Q).transpose();
	
	// loop over all integration points
	vec3d n0e, n0a, nt;
	double In, Wll;
	const double eps = 0;
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	tens4ds N4;
	tens4ds c;
	c.zero();
	
	// right Cauchy-Green tensor: C = Ft*F
	mat3ds C = (F.transpose()*F).sym();
	
	const int nint = (m_nres == 0? NSTL  : NSTH  );
	
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in reference configuration
		n0e.x = m_cth[n]*m_sph[n];
		n0e.y = m_sth[n]*m_sph[n];
		n0e.z = m_cph[n];
		
		// Calculate In = n0e*C*n0e
		In = n0e*(C*n0e);
		
		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// get the global spatial fiber direction in current configuration
			nt = (F*n0e)/sqrt(In);
			
			// calculate the outer product of nt
			N2 = dyad(nt);
			
			// get the local material fiber direction in reference configuration
			n0a = QT*n0e;
			
			// calculate strain energy derivative
			Wll = m_beta*(m_beta-1.0)*m_ksi*pow(In - 1.0, m_beta-2.0);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(4.0/J*In*In*Wll*m_w[n]);
		}
	}
	
	return c;
}
