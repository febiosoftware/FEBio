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
#include "FEEFDUncoupled.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FEEFDUncoupled
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEEFDUncoupled, FEUncoupledMaterial)
	ADD_PARAMETER(m_beta, 3, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETER(m_ksi , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi" )->setUnits(UNIT_PRESSURE);

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEEFDUncoupled::FEEFDUncoupled(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

//-----------------------------------------------------------------------------
mat3ds FEEFDUncoupled::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	// deviatoric deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, Wl;
	const double eps = 0;
	mat3ds s;
	s.zero();

	const int li[4][3] = {
		{ 1, 1, 1},
		{-1, 1, 1},
		{-1,-1, 1},
		{ 1,-1, 1}
	};

	double ksi[3] = { m_ksi[0](mp), m_ksi[1](mp), m_ksi[2](mp) };
	double beta[3] = { m_beta[0](mp), m_beta[1](mp), m_beta[2](mp) };

	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];

		// calculate material coefficients
		// TODO: There is an obvious optimization opportunity here, since the values of ksi
		//       and beta can be precalculated and reused. I have not done this yet since I
		//       need to figure out how to initialize the material parameters for each time
		//       step (instead of once at the start) in case the values depend on load curves.
		double ksi_n  = 1.0 / sqrt(SQR(n0a.x / ksi [0]) + SQR(n0a.y / ksi [1]) + SQR(n0a.z / ksi [2]));
		double beta_n = 1.0 / sqrt(SQR(n0a.x / beta[0]) + SQR(n0a.y / beta[1]) + SQR(n0a.z / beta[2]));

		// loop over the four quadrants
		for (int l=0; l<4; ++l)
		{
			n0q = vec3d(li[l][0]*n0a.x, li[l][1]*n0a.y, li[l][2]*n0a.z);

			// rotate to reference configuration
			n0e = Q*n0q;

			// get the global spatial fiber direction in current configuration
			nt = F*n0e;

			// Calculate In = n0e*C*n0e
			In = nt*nt;
		
			// only take fibers in tension into consideration
			if (In > 1. + eps)
			{
				// calculate the outer product of nt
				mat3ds N = dyad(nt);
			
				// calculate strain energy derivative
				Wl = beta_n *ksi_n *pow(In - 1.0, beta_n -1.0);
			
				// calculate the stress
				s += N*(Wl*wn);
			}
		}
	}
	// don't forget to multiply by two to include the other half-sphere
	return s.dev()*(4.0/J);
}


//-----------------------------------------------------------------------------
tens4ds FEEFDUncoupled::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	// deviatoric deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// loop over all integration points
	vec3d n0e, n0a, nt;
	double In, Wl, Wll;
	const double eps = 0;
	mat3ds s;
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	mat3dd I(1);
	tens4ds N4;
	tens4ds c;

	s.zero();
	c.zero();

	const int li[4][3] = {
		{ 1, 1, 1},
		{-1, 1, 1},
		{-1,-1, 1},
		{ 1,-1, 1}
	};
	
	double ksi[3] = { m_ksi[0](mp), m_ksi[1](mp), m_ksi[2](mp) };
	double beta[3] = { m_beta[0](mp), m_beta[1](mp), m_beta[2](mp) };

	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];

		// calculate material coefficients
		double ksi_n  = 1.0 / sqrt(SQR(n0a.x / ksi [0]) + SQR(n0a.y / ksi [1]) + SQR(n0a.z / ksi [2]));
		double beta_n = 1.0 / sqrt(SQR(n0a.x / beta[0]) + SQR(n0a.y / beta[1]) + SQR(n0a.z / beta[2]));

		for (int l=0; l<4; ++l)
		{
			vec3d n0q = vec3d(li[l][0]*n0a.x, li[l][1]*n0a.y, li[l][2]*n0a.z);

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
				Wl = beta_n *ksi_n *pow(In - 1.0, beta_n -1.0);
				Wll = beta_n *(beta_n -1.0)*ksi_n *pow(In - 1.0, beta_n -2.0);
			
				// calculate the outer product of nt
				N2 = dyad(nt);
				N4 = dyad1s(N2);

				// calculate the stress
				s += N2*(2.0/J*Wl*wn);

				// calculate tangent
				c += N4*(4.0/J*Wll*wn);
			}
		}
	}

	// don't forget to multiply by two to include the other half-sphere
	s *= 2.0;
	c *= 2.0;
	
	// This is the final value of the elasticity tensor
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	c += ((I4+IxI/3.0)*s.tr() - dyad1s(I,s))*(2./3.)
	- (ddots(IxI, c)-IxI*(c.tr()/3.))/3.;
	
	return c;
}

//-----------------------------------------------------------------------------
//! calculate deviatoric strain energy density
double FEEFDUncoupled::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	// deviatoric deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, W;
	const double eps = 0;
    
	const int li[4][3] = {
		{ 1, 1, 1},
		{-1, 1, 1},
		{-1,-1, 1},
		{ 1,-1, 1}
	};
    
	double ksi[3] = { m_ksi[0](mp), m_ksi[1](mp), m_ksi[2](mp) };
	double beta[3] = { m_beta[0](mp), m_beta[1](mp), m_beta[2](mp) };

	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];
        
		// calculate material coefficients
		// TODO: There is an obvious optimization opportunity here, since the values of ksi
		//       and beta can be precalculated and reused. I have not done this yet since I
		//       need to figure out how to initialize the material parameters for each time
		//       step (instead of once at the start) in case the values depend on load curves.
		double ksi_n  = 1.0 / sqrt(SQR(n0a.x / ksi [0]) + SQR(n0a.y / ksi [1]) + SQR(n0a.z / ksi [2]));
		double beta_n = 1.0 / sqrt(SQR(n0a.x / beta[0]) + SQR(n0a.y / beta[1]) + SQR(n0a.z / beta[2]));
        
		// loop over the four quadrants
		for (int l=0; l<4; ++l)
		{
			n0q = vec3d(li[l][0]*n0a.x, li[l][1]*n0a.y, li[l][2]*n0a.z);
            
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
				W = ksi_n *pow(In - 1.0, beta_n);
				sed += W*wn;
			}
		}
	}
    
	// don't forget to multiply by two to include the other half-sphere
    return sed*2.0;
}
