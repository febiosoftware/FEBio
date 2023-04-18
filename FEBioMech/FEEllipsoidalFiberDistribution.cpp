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
#include "FEEllipsoidalFiberDistribution.h"
#include "FEContinuousFiberDistribution.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FEEllipsoidalFiberDistribution
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEEllipsoidalFiberDistribution, FEElasticMaterial)
	ADD_PARAMETER(m_beta, 3, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETER(m_ksi , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi" )->setUnits(UNIT_PRESSURE);

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
mat3ds FEEllipsoidalFiberDistribution::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, Wl;
	const double eps = 0;
	mat3ds s;
	s.zero();

	// calculate material coefficients
	double ksi0 = m_ksi[0](mp);
	double ksi1 = m_ksi[1](mp);
	double ksi2 = m_ksi[2](mp);

	double beta0 = m_beta[0](mp);
	double beta1 = m_beta[1](mp);
	double beta2 = m_beta[2](mp);


	for (int n=0; n<MAX_INT; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];

		// calculate material coefficients
		double ksi  = 1.0 / sqrt(SQR(n0a.x /  ksi0) + SQR(n0a.y /  ksi1) + SQR(n0a.z /  ksi2));
		double beta = 1.0 / sqrt(SQR(n0a.x / beta0) + SQR(n0a.y / beta1) + SQR(n0a.z / beta2));

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
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);

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
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);

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
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);

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
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);

			// calculate the stress
			s += dyad(nt)*(Wl*wn);
		}
	}

	// we multiply by two to add contribution from other half-sphere
	return s*(4.0/J);
}

//-----------------------------------------------------------------------------
tens4ds FEEllipsoidalFiberDistribution::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;

	// get the local coordinate systems
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
	
	// calculate material coefficients
	double ksi0 = m_ksi[0](mp);
	double ksi1 = m_ksi[1](mp);
	double ksi2 = m_ksi[2](mp);

	double beta0 = m_beta[0](mp);
	double beta1 = m_beta[1](mp);
	double beta2 = m_beta[2](mp);

	for (int n=0; n<MAX_INT; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];

		// calculate material coefficients
		double ksi  = 1.0 / sqrt(SQR(n0a.x /  ksi0) + SQR(n0a.y /  ksi1) + SQR(n0a.z /  ksi2));
		double beta = 1.0 / sqrt(SQR(n0a.x / beta0) + SQR(n0a.y / beta1) + SQR(n0a.z / beta2));

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
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);

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
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
			
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
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
			
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
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(Wll*wn);
		}
	}

	// multiply by two to integrate over other half of sphere
	return c*(2.0*4.0/J);
}

//-----------------------------------------------------------------------------
double FEEllipsoidalFiberDistribution::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// loop over all integration points
	vec3d n0e, n0a, n0q, nt;
	double In, W;
	const double eps = 0;

	// calculate material coefficients
	double ksi0 = m_ksi[0](mp);
	double ksi1 = m_ksi[1](mp);
	double ksi2 = m_ksi[2](mp);

	double beta0 = m_beta[0](mp);
	double beta1 = m_beta[1](mp);
	double beta2 = m_beta[2](mp);
    
	const int nint = 45;
	for (int n=0; n<nint; ++n)
	{
		// set the global fiber direction in material coordinate system
		n0a.x = XYZ2[n][0];
		n0a.y = XYZ2[n][1];
		n0a.z = XYZ2[n][2];
		double wn = XYZ2[n][3];
        
		// calculate material coefficients
		double ksi = 1.0 / sqrt(SQR(n0a.x / ksi0) + SQR(n0a.y / ksi1) + SQR(n0a.z / ksi2));
		double beta = 1.0 / sqrt(SQR(n0a.x / beta0) + SQR(n0a.y / beta1) + SQR(n0a.z / beta2));

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
			// calculate strain energy
			W = ksi*pow(In - 1.0, beta);
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
			// calculate strain energy
			W = ksi*pow(In - 1.0, beta);
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
			// calculate strain energy
			W = ksi*pow(In - 1.0, beta);
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
			// calculate strain energy
			W = ksi*pow(In - 1.0, beta);
			sed += W*wn;
		}
	}
    
	// multiply by two to integrate over other half of sphere
    return sed*2.0;
}

//-----------------------------------------------------------------------------
// FEEllipsoidalFiberDistributionOld
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEEllipsoidalFiberDistributionOld, FEElasticMaterial)
	ADD_PARAMETER(m_beta, 3, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETER(m_ksi , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// FEEllipsoidalFiberDistributionOld
//-----------------------------------------------------------------------------

FEEllipsoidalFiberDistributionOld::FEEllipsoidalFiberDistributionOld(FEModel* pfem) : FEElasticMaterial(pfem)
{ 
	m_nres = 0; 
}

//-----------------------------------------------------------------------------
bool FEEllipsoidalFiberDistributionOld::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

	InitIntegrationRule();

	return true;
}

//-----------------------------------------------------------------------------
void FEEllipsoidalFiberDistributionOld::InitIntegrationRule()
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
}

//-----------------------------------------------------------------------------
void FEEllipsoidalFiberDistributionOld::Serialize(DumpStream& ar)
{
	FEElasticMaterial::Serialize(ar);
	if (ar.IsShallow() == false)
	{
		if (ar.IsSaving())
		{
			ar << m_nres;
		}
		else
		{
			ar >> m_nres;
			InitIntegrationRule();
		}
	}
}

//-----------------------------------------------------------------------------
mat3ds FEEllipsoidalFiberDistributionOld::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);
	mat3d QT = Q.transpose();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;

	// loop over all integration points
	double ksi, beta;
	vec3d n0e, n0a, nt;
	double In, Wl;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds s;
	s.zero();

	// material coefficients
	double ksi0 = m_ksi[0](mp);
	double ksi1 = m_ksi[1](mp);
	double ksi2 = m_ksi[2](mp);

	double beta0 = m_beta[0](mp);
	double beta1 = m_beta[1](mp);
	double beta2 = m_beta[2](mp);

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
			
			// calculate material coefficients
			ksi  = 1.0 / sqrt(SQR(n0a.x /  ksi0) + SQR(n0a.y / ksi1 ) + SQR(n0a.z / ksi2 ));
			beta = 1.0 / sqrt(SQR(n0a.x / beta0) + SQR(n0a.y / beta1) + SQR(n0a.z / beta2));
			
			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);
			
			// calculate the stress
			s += N*(2.0/J*In*Wl*m_w[n]);
		}
		
	}
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEEllipsoidalFiberDistributionOld::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);
	mat3d QT = Q.transpose();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;

	// loop over all integration points
	double ksi, beta;
	vec3d n0e, n0a, nt;
	double In, Wll;
	const double eps = 0;
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	tens4ds N4;
	tens4ds c;
	c.zero();
	
	// material coefficients
	double ksi0 = m_ksi[0](mp);
	double ksi1 = m_ksi[1](mp);
	double ksi2 = m_ksi[2](mp);

	double beta0 = m_beta[0](mp);
	double beta1 = m_beta[1](mp);
	double beta2 = m_beta[2](mp);

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
			
			// calculate material coefficients in local fiber direction
			ksi  = 1.0 / sqrt(SQR(n0a.x /  ksi0) + SQR(n0a.y / ksi1 ) + SQR(n0a.z / ksi2 ));
			beta = 1.0 / sqrt(SQR(n0a.x / beta0) + SQR(n0a.y / beta1) + SQR(n0a.z / beta2));
			
			// calculate strain energy derivative
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
			
			N2 = dyad(nt);
			N4 = dyad1s(N2);
			
			c += N4*(4.0/J*In*In*Wll*m_w[n]);
		}
	}
	
	return c;
}

//-----------------------------------------------------------------------------
double FEEllipsoidalFiberDistributionOld::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);
	mat3d QT = Q.transpose();
	
	// deformation gradient
	mat3d &F = pt.m_F;
    
	// loop over all integration points
	double ksi, beta;
	vec3d n0e, n0a, nt;
	double In, W;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
    
	const int nint = (m_nres == 0? NSTL  : NSTH  );
    
	// material coefficients
	double ksi0 = m_ksi[0](mp);
	double ksi1 = m_ksi[1](mp);
	double ksi2 = m_ksi[2](mp);

	double beta0 = m_beta[0](mp);
	double beta1 = m_beta[1](mp);
	double beta2 = m_beta[2](mp);

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
			
			// get the local material fiber direction in reference configuration
			n0a = QT*n0e;
			
			// calculate material coefficients
			ksi  = 1.0 / sqrt(SQR(n0a.x /  ksi0) + SQR(n0a.y / ksi1 ) + SQR(n0a.z / ksi2 ));
			beta = 1.0 / sqrt(SQR(n0a.x / beta0) + SQR(n0a.y / beta1) + SQR(n0a.z / beta2));
			
			// calculate strain energy density
			W = ksi*pow(In - 1.0, beta);
			sed += W*m_w[n];
		}
		
	}
    
    return sed;
}
