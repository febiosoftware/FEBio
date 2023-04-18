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
#include "FEFiberNeoHookean.h"
#include "FEFiberEFDNeoHookean.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberEFDNeoHookean, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
	ADD_PARAMETER(m_a, 3, "a");
	ADD_PARAMETER(m_ac, "active_contraction")->setUnits(UNIT_PRESSURE);

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

#ifndef SQR
	#define SQR(x) ((x)*(x))
#endif

//////////////////////////////////////////////////////////////////////
// FEFiberNeoHookean
//////////////////////////////////////////////////////////////////////


FEFiberEFDNeoHookean::FEFiberEFDNeoHookean(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_nres = 0;
	m_a[0] = m_a[1] = m_a[2] = 0;
	m_ac = 0;
}

bool FEFiberEFDNeoHookean::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

    m_bfirst = true;
    
	if (m_bfirst)
	{
		// select the integration rule
		const int nint    = (m_nres == 0? NSTL  : NSTH  );
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
        
		m_bfirst = false;
	}

	return true;
}

mat3ds FEFiberEFDNeoHookean::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.m_F;
	double detF = pt.m_J;
	double detFi = 1.0/detF;
	double lndetF = log(detF);

	// calculate left Cauchy-Green tensor
	// (we commented out the matrix components we do not need)
	double b[3][3];

	b[0][0] = F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2];
	b[0][1] = F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2];
	b[0][2] = F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2];

//	b[1][0] = F[1][0]*F[0][0]+F[1][1]*F[0][1]+F[1][2]*F[0][2];
	b[1][1] = F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2];
	b[1][2] = F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2];

//	b[2][0] = F[2][0]*F[0][0]+F[2][1]*F[0][1]+F[2][2]*F[0][2];
//	b[2][1] = F[2][0]*F[1][0]+F[2][1]*F[1][1]+F[2][2]*F[1][2];
	b[2][2] = F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2];

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	// calculate stress
	mat3ds s;

	s.xx() = (mu*(b[0][0] - 1) + lam*lndetF)*detFi;
	s.yy() = (mu*(b[1][1] - 1) + lam*lndetF)*detFi;
	s.zz() = (mu*(b[2][2] - 1) + lam*lndetF)*detFi;
	s.xy() = mu*b[0][1]*detFi;
	s.yz() = mu*b[1][2]*detFi;
	s.xz() = mu*b[0][2]*detFi;

	// --- F I B E R   C O N T R I B U T I O N ---
	if (m_ac)
	{
		// select the integration rule
		const int nint    = (m_nres == 0? NSTL  : NSTH  );
		const double* phi = (m_nres == 0? PHIL  : PHIH  );
		const double* the = (m_nres == 0? THETAL: THETAH);
		const double* w   = (m_nres == 0? AREAL : AREAH );

		// get the local coordinate systems
		mat3d Q = GetLocalCS(mp);

		// loop over all integration points
		double nr[3], n0[3], nt[3];
		double at;
		for (int n=0; n<nint; ++n)
		{
			// set the local fiber direction
			nr[0] = m_cth[n]*m_sph[n];
			nr[1] = m_sth[n]*m_sph[n];
			nr[2] = m_cph[n];

			// get the global material fiber direction
			n0[0] = Q[0][0]*nr[0] + Q[0][1]*nr[1] + Q[0][2]*nr[2];
			n0[1] = Q[1][0]*nr[0] + Q[1][1]*nr[1] + Q[1][2]*nr[2];
			n0[2] = Q[2][0]*nr[0] + Q[2][1]*nr[1] + Q[2][2]*nr[2];

			// get the global spatial fiber direction
			nt[0] = F[0][0]*n0[0] + F[0][1]*n0[1] + F[0][2]*n0[2];
			nt[1] = F[1][0]*n0[0] + F[1][1]*n0[1] + F[1][2]*n0[2];
			nt[2] = F[2][0]*n0[0] + F[2][1]*n0[1] + F[2][2]*n0[2];

			// add active contraction stuff
			at = m_ac *sqrt(SQR(m_a[0]*nr[0]) + SQR(m_a[1]*nr[1]) + SQR(m_a[2]*nr[2]));

			s.xx() += at*nt[0]*nt[0];
			s.yy() += at*nt[1]*nt[1];
			s.zz() += at*nt[2]*nt[2];
			s.xy() += at*nt[0]*nt[1];
			s.yz() += at*nt[1]*nt[2];
			s.xz() += at*nt[0]*nt[2];
		}
	}

	return s;
}

tens4ds FEFiberEFDNeoHookean::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.m_F;
	double detF = pt.m_J;

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	double lam1 = lam / detF;
	double mu1  = (mu - lam*log(detF)) / detF;
	
	double D[6][6] = {0};
	D[0][0] = lam1+2.*mu1; D[0][1] = lam1       ; D[0][2] = lam1       ;
	D[1][0] = lam1       ; D[1][1] = lam1+2.*mu1; D[1][2] = lam1       ;
	D[2][0] = lam1       ; D[2][1] = lam1       ; D[2][2] = lam1+2.*mu1;
	D[3][3] = mu1;
	D[4][4] = mu1;
	D[5][5] = mu1;

	return tens4ds(D);
}

double FEFiberEFDNeoHookean::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	double J = pt.m_J;
	double lnJ = log(J);
    
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
    
    double I1 = b.tr();
    
	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);
    
	// calculate strain energy density
	double sed = mu*((I1-3)/2.0 - lnJ) + lam*lnJ*lnJ/2.0;
    
	// --- F I B E R   C O N T R I B U T I O N ---
    // active contraction does not contribute to the strain energy density
    
	return sed;
}
