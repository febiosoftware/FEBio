/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FESphericalFiberDistribution.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// define the material parameters
BEGIN_FECORE_CLASS(FESphericalFiberDistribution, FEElasticMaterial)
	ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
	ADD_PARAMETER(m_beta , FE_RANGE_GREATER_OR_EQUAL(2.0), "beta" );
	ADD_PARAMETER(m_ksi  , FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi"  );

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// FESphericalFiberDistribution
//-----------------------------------------------------------------------------

FESphericalFiberDistribution::FESphericalFiberDistribution(FEModel* pfem) : FEElasticMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
mat3ds FESphericalFiberDistribution::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// get local coordinates
	mat3d Q = GetLocalCS(mp);

	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, Wl;
	const double eps = 0;
	mat3ds s;
	s.zero();
	
	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];
		
		// calculate material coefficients
		double ksi  = m_ksi(mp);
		double alpha = m_alpha(mp);
		double beta = m_beta(mp);
		
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
tens4ds FESphericalFiberDistribution::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// get local coordinates
	mat3d Q = GetLocalCS(mp);

	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, Wll;
	const double eps = 0;
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	tens4ds N4;
	tens4ds c;
	c.zero();
	
	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];
		
		// calculate material coefficients
		double ksi  = m_ksi(mp);
		double alpha = m_alpha(mp);
		double beta = m_beta(mp);
		
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
double FESphericalFiberDistribution::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// get local coordinates
	mat3d Q = GetLocalCS(mp);

	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, W;
	const double eps = 0;
	double sed = 0.0;
	
	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];
		
		// calculate material coefficients
		double ksi  = m_ksi(mp);
		double alpha = m_alpha(mp);
		double beta = m_beta(mp);
		
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
            if (alpha > 0)
                W = ksi/alpha*(exp(alpha*pow(In - 1.0, beta))-1);
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
            if (alpha > 0)
                W = ksi/alpha*(exp(alpha*pow(In - 1.0, beta))-1);
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
            if (alpha > 0)
                W = ksi/alpha*(exp(alpha*pow(In - 1.0, beta))-1);
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
            if (alpha > 0)
                W = ksi/alpha*(exp(alpha*pow(In - 1.0, beta))-1);
            else
                W = ksi*pow(In - 1.0, beta);
			
			// add to total sed
			sed += W*wn;
		}
	}
	
	// we multiply by two to add contribution from other half-sphere
	return sed*2.0;
}
