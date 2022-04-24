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
#include "FE2DFiberNeoHookean.h"

// define the material parameters
BEGIN_FECORE_CLASS(FE2DFiberNeoHookean, FEElasticMaterial)
	ADD_PARAMETER(m_E, "E")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_v, "v");
	ADD_PARAMETER(m_a, 2, "a");
	ADD_PARAMETER(m_ac, "active_contraction")->setUnits(UNIT_PRESSURE);
	
	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);
END_FECORE_CLASS();

double FE2DFiberNeoHookean::m_cth[FE2DFiberNeoHookean::NSTEPS];
double FE2DFiberNeoHookean::m_sth[FE2DFiberNeoHookean::NSTEPS];

//////////////////////////////////////////////////////////////////////
// FE2DFiberNeoHookean
//////////////////////////////////////////////////////////////////////

FE2DFiberNeoHookean::FE2DFiberNeoHookean(FEModel* pfem) : FEElasticMaterial(pfem)
{
	static bool bfirst = true;

	if (bfirst)
	{
		double ph;
		for (int n=0; n<NSTEPS; ++n)
		{
			ph = 2.0*PI*n / (double) NSTEPS;
			m_cth[n] = cos(ph);
			m_sth[n] = sin(ph);
		}
		bfirst = false;
	}

	m_ac = 0;
	m_a[0] = m_a[1] = 0;
}

mat3ds FE2DFiberNeoHookean::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.m_F;
	double detF = pt.m_J;
	double detFi = 1.0/detF;
	double lndetF = log(detF);

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

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

	// material parameters
	double E = m_E(mp);
	double v = m_v(mp);

	// lame parameters
	double lam = v*E/((1+v)*(1-2*v));
	double mu  = 0.5*E/(1+v);

	// calculate stress
	mat3ds s;

	// --- M A T R I X   C O N T R I B U T I O N ---
	s.xx() = (mu*(b[0][0] - 1) + lam*lndetF)*detFi;
	s.yy() = (mu*(b[1][1] - 1) + lam*lndetF)*detFi;
	s.zz() = (mu*(b[2][2] - 1) + lam*lndetF)*detFi;
	s.xy() = mu*b[0][1]*detFi;
	s.yz() = mu*b[1][2]*detFi;
	s.xz() = mu*b[0][2]*detFi;

	// --- F I B E R   C O N T R I B U T I O N ---
	// NOTE: we have only implemented the active contraction model for this material
	// There is no passive fiber stress.
	if (m_ac > 0)
	{
		double wa = 1.0 / (double) NSTEPS;
		vec3d a0, a, v;
		double at;
		for (int n=0; n<NSTEPS; ++n)
		{
			// calculate the local material fiber vector
			v.x = m_cth[n];
			v.y = m_sth[n];
			v.z = 0;

			// calculate the global material fiber vector
			a0 = Q*v;

			// calculate the global spatial fiber vector
			a.x = F[0][0]*a0.x + F[0][1]*a0.y + F[0][2]*a0.z;
			a.y = F[1][0]*a0.x + F[1][1]*a0.y + F[1][2]*a0.z;
			a.z = F[2][0]*a0.x + F[2][1]*a0.y + F[2][2]*a0.z;

			// normalize material axis and store fiber stretch
			a.unit();

			// add active contraction stuff
			at = wa*m_ac *sqrt((m_a[0]*v.x)*(m_a[0]*v.x) + (m_a[1]*v.y)*(m_a[1]*v.y));

			s.xx() += at*a.x*a.x;
			s.yy() += at*a.y*a.y;
			s.zz() += at*a.z*a.z;
			s.xy() += at*a.x*a.y;
			s.yz() += at*a.y*a.z;
			s.xz() += at*a.x*a.z;
		}
	}

	return s;
}

tens4ds FE2DFiberNeoHookean::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.m_F;
	double detF = pt.m_J;

	// material parameters
	double E = m_E(mp);
	double v = m_v(mp);

	// lame parameters
	double lam = v*E/((1+v)*(1-2*v));
	double mu  = 0.5*E/(1+v);

	double lam1 = lam / detF;
	double mu1  = (mu - lam*log(detF)) / detF;

	// --- M A T R I X   C O N T R I B U T I O N ---
	double D[6][6] = {0};
	D[0][0] = lam1+2.*mu1; D[0][1] = lam1       ; D[0][2] = lam1       ;
	D[1][0] = lam1       ; D[1][1] = lam1+2.*mu1; D[1][2] = lam1       ;
	D[2][0] = lam1       ; D[2][1] = lam1       ; D[2][2] = lam1+2.*mu1;
	D[3][3] = mu1;
	D[4][4] = mu1;
	D[5][5] = mu1;

	// --- F I B E R   C O N T R I B U T I O N ---
	// NOTE: I commented the fiber stiffness out since I think it will lead to
	// a nonsymmetric D matrix and I can't deal with that yet. Besides, most
	// problems seem to be running just fine without this contribution.
/*
	if (m_ac)
	{
		// Next, we add the fiber contribution. Since the fibers lie
		// randomly perpendicular to the transverse axis, we need
		// to integrate over that plane
		const double PI = 4.0*atan(1.0);
		double lam, at, In;
		vec3d a0, a, v;
		double wa = 1.0 / (double) NSTEPS;
		for (int n=0; n<NSTEPS; ++n)
		{
			// calculate the local material fiber vector
			v.x = m_cth[n];
			v.y = m_sth[n];
			v.z = 0;

			// calculate the global material fiber vector
			a0 = pt.Q*v;

			// calculate the global spatial fiber vector
			a.x = F[0][0]*a0.x + F[0][1]*a0.y + F[0][2]*a0.z;
			a.y = F[1][0]*a0.x + F[1][1]*a0.y + F[1][2]*a0.z;
			a.z = F[2][0]*a0.x + F[2][1]*a0.y + F[2][2]*a0.z;

			// normalize material axis and store fiber stretch
			lam = a.unit();

			// add active contraction stuff
			at = wa*m_ac *sqrt((m_a[0]*v.x)*(m_a[0]*v.x) + (m_a[1]*v.y)*(m_a[1]*v.y));
			In = lam*lam;

			D[0][0] += at*(a.x*a.x - 2.0*a.x*a.x*a.x*a.x);	// c(0,0,0,0)
			D[0][1] += at*(a.x*a.x - 2.0*a.x*a.x*a.y*a.y);	// c(0,0,1,1)
			D[0][2] += at*(a.x*a.x - 2.0*a.x*a.x*a.z*a.z); // c(0,0,2,2)
			D[0][3] -= at*(2.0*a.x*a.x*a.x*a.y); // c(0,0,0,1)
			D[0][4] -= at*(2.0*a.x*a.x*a.y*a.z); // c(0,0,1,2)
			D[0][5] -= at*(2.0*a.x*a.x*a.x*a.z); // c(0,0,0,2)

			D[1][1] += at*(a.y*a.y - 2.0*a.y*a.y*a.y*a.y); // c(1,1,1,1)
			D[1][2] += at*(a.y*a.y - 2.0*a.y*a.y*a.z*a.z); // c(1,1,2,2)
			D[1][3] -= at*(2.0*a.y*a.y*a.x*a.y); // c(1,1,0,1)
			D[1][4] -= at*(2.0*a.y*a.y*a.y*a.z); // c(1,1,1,2)
			D[1][5] -= at*(2.0*a.y*a.y*a.x*a.z); // c(1,1,0,2)

			D[2][2] += at*(a.z*a.z - 2.0*a.z*a.z*a.z*a.z); // c(2,2,2,2)
			D[2][3] -= at*(2.0*a.z*a.z*a.x*a.y); // c(2,2,0,1)
			D[2][4] -= at*(2.0*a.z*a.z*a.y*a.z); // c(2,2,1,2)
			D[2][5] -= at*(2.0*a.z*a.z*a.x*a.z); // c(2,2,0,2)

			D[3][3] -= at*(2.0*a.x*a.y*a.x*a.y); // c(0,1,0,1)
			D[3][4] -= at*(2.0*a.x*a.y*a.y*a.z); // c(0,1,1,2)
			D[3][5] -= at*(2.0*a.x*a.y*a.x*a.z); // c(0,1,0,2)

			D[4][4] -= at*(2.0*a.y*a.z*a.y*a.z); // c(1,2,1,2)
			D[4][5] -= at*(2.0*a.y*a.z*a.x*a.z); // c(1,2,0,2)

			D[5][5] -= at*(2.0*a.x*a.z*a.x*a.z); // c(0,2,0,2)
		}
	}
*/
	return tens4ds(D);
}
