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
#include <limits>
#include "FETCNonlinearOrthotropic.h"

// define the material parameters
BEGIN_FECORE_CLASS(FETCNonlinearOrthotropic, FEUncoupledMaterial)
	ADD_PARAMETER(m_c1, "c1");
	ADD_PARAMETER(m_c2, "c2");
	ADD_PARAMETER(m_beta, 3, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETER(m_ksi , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi" );

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);
END_FECORE_CLASS();

//////////////////////////////////////////////////////////////////////
// FETCNonlinearOrthotropic
//////////////////////////////////////////////////////////////////////

FETCNonlinearOrthotropic::FETCNonlinearOrthotropic(FEModel* pfem) : FEUncoupledMaterial(pfem) 
{
	m_c1 = 0.0;
	m_c2 = 0.0;
	m_beta[0] = m_beta[1] = m_beta[2] = 0.0;
	m_ksi[0] = m_ksi[1] = m_ksi[2] = 0.0;

	m_epsf = 0.0;
}

//-----------------------------------------------------------------------------
// Calculate the deviatoric Cauchy stress
mat3ds FETCNonlinearOrthotropic::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Ji = 1.0 / J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double twoJi = 2.0*Ji;

	// invariants of B
	double I1, I2, I4a, I4b, I4c;

	// strain energy derivatives
	double W1, W2, W4a, W4b, W4c;

	// current local material axis
	vec3d a0, b0, c0, a, b, c;
	double la, lb, lc, lat, lbt, lct;

	double w1pw2i1; // = W1 + W2*I1

	const double third = 1.0/3.0;

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber directions
	a0.x = Q[0][0]; b0.x = Q[0][1]; c0.x = Q[0][2];
	a0.y = Q[1][0]; b0.y = Q[1][1]; c0.y = Q[1][2];
	a0.z = Q[2][0]; b0.z = Q[2][1]; c0.z = Q[2][2];

	// calculate the current material axes lam*a = F*a0;
	a = F*a0;
	b = F*b0;
	c = F*c0;

	// normalize material axis and store fiber stretch
	la  = a.unit();
	lat = la*Jm13; // i.e. lambda tilde

	lb  = b.unit();
	lbt = lb*Jm13;

	lc  = c.unit();
	lct = lc*Jm13;

	// get deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	I1 = B.tr();
	I2 = 0.5*(I1*I1 - B2.tr());
	I4a = lat*lat;
	I4b = lbt*lbt;
	I4c = lct*lct;

	// --- put strain energy derivatives here ---
	// Wi = dW/dIi
	W1 = m_c1;
	W2 = m_c2;

	// fiber a
	if (lat > 1)
	{
		double lati = 1.0/lat;
		double Wl;
		Wl  = m_beta[0]*m_ksi[0]*pow((lat - 1.0), m_beta[0]-1.0);
		W4a = 0.5*lati*Wl;
	}
	else 
	{
		W4a = 0;
	}

	// fiber b
	if (lbt > 1)
	{
		double lbti = 1.0/lbt;
		double Wl;
		Wl  = m_beta[1]*m_ksi[1]*pow((lbt - 1.0), m_beta[1]-1.0);
		W4b = 0.5*lbti*Wl;
	}
	else 
	{
		W4b = 0;
	}

	// fiber c
	if (lct > 1)
	{
		double lcti = 1.0/lct;
		double Wl;
		Wl  = m_beta[2]*m_ksi[2]*pow((lct - 1.0), m_beta[2]-1.0);
		W4c = 0.5*lcti*Wl;
	}
	else 
	{
		W4c = 0;
	}
	// ---

	// calculate T = F*dW/dC*Ft
	// (we commented out the matrix components we do not need)
	w1pw2i1 = W1 + W2*I1;
	
	mat3ds AxA = dyad(a);
	mat3ds BxB = dyad(b);
	mat3ds CxC = dyad(c);

	mat3ds T = B*w1pw2i1 - B2*W2 + AxA*(W4a*I4a) + BxB*(W4b*I4b) + CxC*(W4c*I4c);

	return T.dev()*(2.0/J);
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric tangent
tens4ds FETCNonlinearOrthotropic::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	const double third = 1.0 / 3.0;

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double Ji = 1.0/J;

	// deviatoric cauchy-stress, trs = trace[s]/3
	double s[3][3], trs;
	mat3ds& es = pt.m_s;
	s[0][0] = es.xx();
	s[1][1] = es.yy();
	s[2][2] = es.zz();
	s[0][1] = s[1][0] = es.xy();
	s[1][2] = s[2][1] = es.yz();
	s[0][2] = s[2][0] = es.xz();

	trs = (s[0][0] + s[1][1] + s[2][2])*third;
	s[0][0] -= trs;	
	s[1][1] -= trs;	
	s[2][2] -= trs;

	// current material axis lam*a = F*a0;
	vec3d a, b, c, a0, b0, c0;
	double la, lb, lc, lat, lbt, lct;

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber directions
	a0.x = Q[0][0]; b0.x = Q[0][1]; c0.x = Q[0][2];
	a0.y = Q[1][0]; b0.y = Q[1][1]; c0.y = Q[1][2];
	a0.z = Q[2][0]; b0.z = Q[2][1]; c0.z = Q[2][2];

	// calculate the current material axes lam*a = F*a0;
	a.x = F[0][0]*a0.x + F[0][1]*a0.y + F[0][2]*a0.z;
	a.y = F[1][0]*a0.x + F[1][1]*a0.y + F[1][2]*a0.z;
	a.z = F[2][0]*a0.x + F[2][1]*a0.y + F[2][2]*a0.z;

	b.x = F[0][0]*b0.x + F[0][1]*b0.y + F[0][2]*b0.z;
	b.y = F[1][0]*b0.x + F[1][1]*b0.y + F[1][2]*b0.z;
	b.z = F[2][0]*b0.x + F[2][1]*b0.y + F[2][2]*b0.z;

	c.x = F[0][0]*c0.x + F[0][1]*c0.y + F[0][2]*c0.z;
	c.y = F[1][0]*c0.x + F[1][1]*c0.y + F[1][2]*c0.z;
	c.z = F[2][0]*c0.x + F[2][1]*c0.y + F[2][2]*c0.z;

	la = a.unit();
	lb = b.unit();
	lc = c.unit();

	// deviatoric stretch
	lat = la*Jm13;
	lbt = lb*Jm13;
	lct = lc*Jm13;

	// deviatoric right Cauchy-Green tensor: C = Ft*F
	double C[3][3];
	C[0][0] = Jm23*(F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]);
	C[0][1] = Jm23*(F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]);
	C[0][2] = Jm23*(F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]);

	C[1][0] = Jm23*(F[0][1]*F[0][0]+F[1][1]*F[1][0]+F[2][1]*F[2][0]);
	C[1][1] = Jm23*(F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]);
	C[1][2] = Jm23*(F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]);

	C[2][0] = Jm23*(F[0][2]*F[0][0]+F[1][2]*F[1][0]+F[2][2]*F[2][0]);
	C[2][1] = Jm23*(F[0][2]*F[0][1]+F[1][2]*F[1][1]+F[2][2]*F[2][1]);
	C[2][2] = Jm23*(F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]);

	// square of C
	// (we commented out the components we don't need)
	double C2[3][3];
	C2[0][0] = C[0][0]*C[0][0]+C[0][1]*C[1][0]+C[0][2]*C[2][0];
//	C2[0][1] = C[0][0]*C[0][1]+C[0][1]*C[1][1]+C[0][2]*C[2][1];
//	C2[0][2] = C[0][0]*C[0][2]+C[0][1]*C[1][2]+C[0][2]*C[2][2];

//	C2[1][0] = C[1][0]*C[0][0]+C[1][1]*C[1][0]+C[1][2]*C[2][0];
	C2[1][1] = C[1][0]*C[0][1]+C[1][1]*C[1][1]+C[1][2]*C[2][1];
//	C2[1][2] = C[1][0]*C[0][2]+C[1][1]*C[1][2]+C[1][2]*C[2][2];

//	C2[2][0] = C[2][0]*C[0][0]+C[2][1]*C[1][0]+C[2][2]*C[2][0];
//	C2[2][1] = C[2][0]*C[0][1]+C[2][1]*C[1][1]+C[2][2]*C[2][1];
	C2[2][2] = C[2][0]*C[0][2]+C[2][1]*C[1][2]+C[2][2]*C[2][2];

	// Invariants of C
	double I1 = C[0][0] + C[1][1] + C[2][2];
	double I2 = 0.5*(I1*I1 - (C2[0][0] + C2[1][1] + C2[2][2]));
	double I4a = lat*lat;
	double I4b = lbt*lbt;
	double I4c = lct*lct;

	// calculate left Cauchy-Green tensor: B = F*Ft
	double B[3][3];
	B[0][0] = Jm23*(F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2]);
	B[0][1] = Jm23*(F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2]);
	B[0][2] = Jm23*(F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2]);

	B[1][0] = Jm23*(F[1][0]*F[0][0]+F[1][1]*F[0][1]+F[1][2]*F[0][2]);
	B[1][1] = Jm23*(F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2]);
	B[1][2] = Jm23*(F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2]);

	B[2][0] = Jm23*(F[2][0]*F[0][0]+F[2][1]*F[0][1]+F[2][2]*F[0][2]);
	B[2][1] = Jm23*(F[2][0]*F[1][0]+F[2][1]*F[1][1]+F[2][2]*F[1][2]);
	B[2][2] = Jm23*(F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2]);

	// calculate square of B
	// (we commented out the components we don't need)
	double B2[3][3];
	B2[0][0] = B[0][0]*B[0][0]+B[0][1]*B[1][0]+B[0][2]*B[2][0];
	B2[0][1] = B[0][0]*B[0][1]+B[0][1]*B[1][1]+B[0][2]*B[2][1];
	B2[0][2] = B[0][0]*B[0][2]+B[0][1]*B[1][2]+B[0][2]*B[2][2];

//	B2[1][0] = B[1][0]*B[0][0]+B[1][1]*B[1][0]+B[1][2]*B[2][0];
	B2[1][1] = B[1][0]*B[0][1]+B[1][1]*B[1][1]+B[1][2]*B[2][1];
	B2[1][2] = B[1][0]*B[0][2]+B[1][1]*B[1][2]+B[1][2]*B[2][2];

//	B2[2][0] = B[2][0]*B[0][0]+B[2][1]*B[1][0]+B[2][2]*B[2][0];
//	B2[2][1] = B[2][0]*B[0][1]+B[2][1]*B[1][1]+B[2][2]*B[2][1];
	B2[2][2] = B[2][0]*B[0][2]+B[2][1]*B[1][2]+B[2][2]*B[2][2];

	// --- put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1, W2, W4a, W4b, W4c, W44a, W44b, W44c;
	W1 = m_c1;
	W2 = m_c2;

	const double eps = m_epsf*std::numeric_limits<double>::epsilon();

	// fiber a
	if (lat >= 1 + eps)
	{
		double lati = 1.0/lat;

		double Wl, Wll;
		Wl  = m_beta[0]*m_ksi[0]*pow((lat - 1.0), m_beta[0]-1.0);
		Wll = (m_beta[0]-1.0)*m_beta[0]*m_ksi[0]*pow((lat - 1.0), m_beta[0]-2.0);
		W4a  = 0.5*lati*Wl;
		W44a = 0.25*lati*lati*(Wll - lati*Wl);
	}
	else 
	{
		W4a  = 0;
		W44a = 0;
	}

	// fiber b
	if (lbt >= 1 + eps)
	{
		double lbti = 1.0/lbt;

		double Wl, Wll;
		Wl  = m_beta[1]*m_ksi[1]*pow((lbt - 1.0), m_beta[1]-1.0);
		Wll = (m_beta[1]-1.0)*m_beta[1]*m_ksi[1]*pow((lbt - 1.0), m_beta[1]-2.0);
		W4b  = 0.5*lbti*Wl;
		W44b = 0.25*lbti*lbti*(Wll - lbti*Wl);
	}
	else 
	{
		W4b  = 0;
		W44b = 0;
	}

	// fiber c
	if (lct >= 1 + eps)
	{
		double lcti = 1.0/lct;

		double Wl, Wll;
		Wl  = m_beta[2]*m_ksi[2]*pow((lct - 1.0), m_beta[2]-1.0);
		Wll = (m_beta[2]-1.0)*m_beta[2]*m_ksi[2]*pow((lct - 1.0), m_beta[2]-2.0);
		W4c  = 0.5*lcti*Wl;
		W44c = 0.25*lcti*lcti*(Wll - lcti*Wl);
	}
	else 
	{
		W4c  = 0;
		W44c = 0;
	}
	// ------------------------------------

	double D[6][6] = {0};

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2 + W4a*I4a + W4b*I4b + W4c*I4c;

	// calculate C:d2WdCdC:C
	double CWWC = 2*I2*W2 + W44a*I4a*I4a + W44b*I4b*I4b + W44c*I4c*I4c;

	// D[0][0] = c(0,0,0,0)
	D[0][0] =  -(4.0/3.0)*s[0][0] + (8.0/9.0)*Ji*WC;
	D[0][0] += 4.0*Ji*W44a*I4a*I4a*a.x*a.x*a.x*a.x;
	D[0][0] += 4.0*Ji*W44b*I4b*I4b*b.x*b.x*b.x*b.x;
	D[0][0] += 4.0*Ji*W44c*I4c*I4c*c.x*c.x*c.x*c.x;
	D[0][0] += (4.0/9.0)*Ji*CWWC;
	D[0][0] -= (8.0/3.0)*Ji*(W2*I1*B[0][0] - W2*B2[0][0] + W44a*I4a*I4a*a.x*a.x + W44b*I4b*I4b*b.x*b.x + W44c*I4c*I4c*c.x*c.x);

	// D[1][1] = c(1,1,1,1)
	D[1][1] = -(4.0/3.0)*s[1][1] + 8.0/9.0*Ji*WC;
	D[1][1] += 4.0*Ji*W44a*I4a*I4a*a.y*a.y*a.y*a.y;
	D[1][1] += 4.0*Ji*W44b*I4b*I4b*b.y*b.y*b.y*b.y;
	D[1][1] += 4.0*Ji*W44c*I4c*I4c*c.y*c.y*c.y*c.y;
	D[1][1] += (4.0/9.0)*Ji*CWWC;
	D[1][1] -= (8.0/3.0)*Ji*(W2*I1*B[1][1] - W2*B2[1][1] + W44a*I4a*I4a*a.y*a.y + W44b*I4b*I4b*b.y*b.y + W44c*I4c*I4c*c.y*c.y);

	// D[2][2] = c(2,2,2,2)
	D[2][2] = -(4.0/3.0)*s[2][2] + (8.0/9.0)*Ji*WC;
	D[2][2] += 4.0*Ji*W44a*I4a*I4a*a.z*a.z*a.z*a.z;
	D[2][2] += 4.0*Ji*W44b*I4b*I4b*b.z*b.z*b.z*b.z;
	D[2][2] += 4.0*Ji*W44c*I4c*I4c*c.z*c.z*c.z*c.z;
	D[2][2] += (4.0/9.0)*Ji*CWWC;
	D[2][2] -= (8.0/3.0)*Ji*(W2*I1*B[2][2] - W2*B2[2][2] + W44a*I4a*I4a*a.z*a.z + W44b*I4b*I4b*b.z*b.z + W44c*I4c*I4c*c.z*c.z);



	// D[0][1] = D[1][0] = c(0,0,1,1)
	D[0][1] = -(2.0/3.0)*(s[0][0] + s[1][1]) - (4.0/9.0)*Ji*WC;
	D[0][1] += 4.0*Ji*W44a*I4a*I4a*a.x*a.x*a.y*a.y;
	D[0][1] += 4.0*Ji*W44b*I4b*I4b*b.x*b.x*b.y*b.y;
	D[0][1] += 4.0*Ji*W44c*I4c*I4c*c.x*c.x*c.y*c.y;
	D[0][1] += 4.0*Ji*W2*(B[0][0]*B[1][1] - B[0][1]*B[0][1]);
	D[0][1] += (4.0/9.0)*Ji*CWWC;
	D[0][1] -= (4.0/3.0)*Ji*(W2*(I1*B[0][0] - B2[0][0]) + W44a*I4a*I4a*a.x*a.x + W44b*I4b*I4b*b.x*b.x + W44c*I4c*I4c*c.x*c.x);
	D[0][1] -= (4.0/3.0)*Ji*(W2*(I1*B[1][1] - B2[1][1]) + W44a*I4a*I4a*a.y*a.y + W44b*I4b*I4b*b.y*b.y + W44c*I4c*I4c*c.y*c.y);

	// D[1][2] = D[2][1] = c(1,1,2,2)
	D[1][2] = -(2.0/3.0)*(s[1][1] + s[2][2]) - (4.0/9.0)*Ji*WC;
	D[1][2] += 4.0*Ji*W44a*I4a*I4a*a.y*a.y*a.z*a.z;
	D[1][2] += 4.0*Ji*W44b*I4b*I4b*b.y*b.y*b.z*b.z;
	D[1][2] += 4.0*Ji*W44c*I4c*I4c*c.y*c.y*c.z*c.z;
	D[1][2] += 4.0*Ji*W2*(B[1][1]*B[2][2] - B[1][2]*B[1][2]);
	D[1][2] += (4.0/9.0)*Ji*CWWC;
	D[1][2] -= (4.0/3.0)*Ji*(W2*(I1*B[1][1] - B2[1][1]) + W44a*I4a*I4a*a.y*a.y + W44b*I4b*I4b*b.y*b.y + W44c*I4c*I4c*c.y*c.y);
	D[1][2] -= (4.0/3.0)*Ji*(W2*(I1*B[2][2] - B2[2][2]) + W44a*I4a*I4a*a.z*a.z + W44b*I4b*I4b*b.z*b.z + W44c*I4c*I4c*c.z*c.z);

	// D[0][2] = D[2][0] = c(0,0,2,2)
	D[0][2] = -(2.0/3.0)*(s[0][0] + s[2][2]) - (4.0/9.0)*Ji*WC;
	D[0][2] += 4.0*Ji*W44a*I4a*I4a*a.x*a.x*a.z*a.z;
	D[0][2] += 4.0*Ji*W44b*I4b*I4b*b.x*b.x*b.z*b.z;
	D[0][2] += 4.0*Ji*W44c*I4c*I4c*c.x*c.x*c.z*c.z;
	D[0][2] += 4.0*Ji*W2*(B[0][0]*B[2][2] - B[0][2]*B[0][2]);
	D[0][2] += (4.0/9.0)*Ji*CWWC;
	D[0][2] -= (4.0/3.0)*Ji*(W2*(I1*B[0][0] - B2[0][0]) + W44a*I4a*I4a*a.x*a.x + W44b*I4b*I4b*b.x*b.x + W44c*I4c*I4c*c.x*c.x);
	D[0][2] -= (4.0/3.0)*Ji*(W2*(I1*B[2][2] - B2[2][2]) + W44a*I4a*I4a*a.z*a.z + W44b*I4b*I4b*b.z*b.z + W44c*I4c*I4c*c.z*c.z);



	// D[3][3] = 0.5*(c(0,1,0,1) + c(0,1,1,0))
	D[3][3] = (2.0/3.0)*Ji*WC;
	D[3][3] += 4.0*Ji*W44a*I4a*I4a*a.x*a.y*a.x*a.y;
	D[3][3] += 4.0*Ji*W44b*I4b*I4b*b.x*b.y*b.x*b.y;
	D[3][3] += 4.0*Ji*W44c*I4c*I4c*c.x*c.y*c.x*c.y;
	D[3][3] += 2.0*Ji*W2*(B[0][1]*B[0][1] - B[0][0]*B[1][1]);

	// D[4][4] = 0.5*(c(1,2,1,2) + c(1,2,2,1))
	D[4][4] = (2.0/3.0)*Ji*WC;
	D[4][4] += 4.0*Ji*W44a*I4a*I4a*a.y*a.z*a.y*a.z;
	D[4][4] += 4.0*Ji*W44b*I4b*I4b*b.y*b.z*b.y*b.z;
	D[4][4] += 4.0*Ji*W44c*I4c*I4c*c.y*c.z*c.y*c.z;
	D[4][4] += 2.0*Ji*W2*(B[1][2]*B[1][2] - B[1][1]*B[2][2]);

	// D[5][5] = 0.5*(c(0,2,0,2) + c(0,2,2,0))
	D[5][5] = (2.0/3.0)*Ji*WC;
	D[5][5] += 4.0*Ji*W44a*I4a*I4a*a.x*a.z*a.x*a.z;
	D[5][5] += 4.0*Ji*W44b*I4b*I4b*b.x*b.z*b.x*b.z;
	D[5][5] += 4.0*Ji*W44c*I4c*I4c*c.x*c.z*c.x*c.z;
	D[5][5] += 2.0*Ji*W2*(B[0][2]*B[0][2] - B[0][0]*B[2][2]);



	// D[0][3] = 0.5*(c(0,0,0,1) + c(0,0,1,0))
	D[0][3] =  -(2.0/3.0)*s[0][1];
	D[0][3] += 4.0*Ji*W44a*I4a*I4a*a.x*a.x*a.x*a.y;
	D[0][3] += 4.0*Ji*W44b*I4b*I4b*b.x*b.x*b.x*b.y;
	D[0][3] += 4.0*Ji*W44c*I4c*I4c*c.x*c.x*c.x*c.y;
	D[0][3] -= (4.0/3.0)*Ji*(W2*(I1*B[0][1] - B2[0][1]) + W44a*I4a*I4a*a.x*a.y + W44b*I4b*I4b*b.x*b.y + W44c*I4c*I4c*c.x*c.y);

	// D[0][4] = 0.5*(c(0,0,1,2) + c(0,0,2,1))
	D[0][4] =  -(2.0/3.0)*s[1][2];
	D[0][4] += 4.0*Ji*W44a*I4a*I4a*a.x*a.x*a.y*a.z;
	D[0][4] += 4.0*Ji*W44b*I4b*I4b*b.x*b.x*b.y*b.z;
	D[0][4] += 4.0*Ji*W44c*I4c*I4c*c.x*c.x*c.y*c.z;
	D[0][4] += 4.0*Ji*W2*(B[0][0]*B[1][2] - B[0][1]*B[0][2]);
	D[0][4] -= (4.0/3.0)*Ji*(W2*(I1*B[1][2] - B2[1][2]) + W44a*I4a*I4a*a.y*a.z + W44b*I4b*I4b*b.y*b.z + W44c*I4c*I4c*c.y*c.z);

	// D[0][5] = 0.5*(c(0,0,0,2) + c(0,0,2,0))
	D[0][5] =  -(2.0/3.0)*s[0][2];
	D[0][5] += 4.0*Ji*W44a*I4a*I4a*a.x*a.x*a.x*a.z;
	D[0][5] += 4.0*Ji*W44b*I4b*I4b*b.x*b.x*b.x*b.z;
	D[0][5] += 4.0*Ji*W44c*I4c*I4c*c.x*c.x*c.x*c.z;
	D[0][5] -= (4.0/3.0)*Ji*(W2*(I1*B[0][2] - B2[0][2]) + W44a*I4a*I4a*a.x*a.z + W44b*I4b*I4b*b.x*b.z + W44c*I4c*I4c*c.x*c.z);

	// D[1][3] = 0.5*(c(1,1,0,1) + c(1,1,1,0))
	D[1][3] =  -(2.0/3.0)*s[0][1];
	D[1][3] += 4.0*Ji*W44a*I4a*I4a*a.y*a.y*a.x*a.y;
	D[1][3] += 4.0*Ji*W44b*I4b*I4b*b.y*b.y*b.x*b.y;
	D[1][3] += 4.0*Ji*W44c*I4c*I4c*c.y*c.y*c.x*c.y;
	D[1][3] -= (4.0/3.0)*Ji*(W2*(I1*B[0][1] - B2[0][1]) + W44a*I4a*I4a*a.x*a.y + W44b*I4b*I4b*b.x*b.y + W44c*I4c*I4c*c.x*c.y);

	// D[1][4] = 0.5*(c(1,1,1,2) + c(1,1,2,1))
	D[1][4] =  -(2.0/3.0)*s[1][2];
	D[1][4] += 4.0*Ji*W44a*I4a*I4a*a.y*a.y*a.y*a.z;
	D[1][4] += 4.0*Ji*W44b*I4b*I4b*b.y*b.y*b.y*b.z;
	D[1][4] += 4.0*Ji*W44c*I4c*I4c*c.y*c.y*c.y*c.z;
	D[1][4] -= (4.0/3.0)*Ji*(W2*(I1*B[1][2] - B2[1][2]) + W44a*I4a*I4a*a.y*a.z + W44b*I4b*I4b*b.y*b.z + W44c*I4c*I4c*c.y*c.z);

	// D[1][5] = 0.5*(c(1,1,0,2) + c(1,1,2,0))
	D[1][5] =  -(2.0/3.0)*s[0][2];
	D[1][5] += 4.0*Ji*W44a*I4a*I4a*a.y*a.y*a.x*a.z;
	D[1][5] += 4.0*Ji*W44b*I4b*I4b*b.y*b.y*b.x*b.z;
	D[1][5] += 4.0*Ji*W44c*I4c*I4c*c.y*c.y*c.x*c.z;
	D[1][5] += 4.0*Ji*W2*(B[1][1]*B[0][2] - B[0][1]*B[1][2]);
	D[1][5] -= (4.0/3.0)*Ji*(W2*(I1*B[0][2] - B2[0][2]) + W44a*I4a*I4a*a.x*a.z + W44b*I4b*I4b*b.x*b.z + W44c*I4c*I4c*c.x*c.z);

	// D[2][3] = 0.5*(c(2,2,0,1) + c(2,2,1,0))
	D[2][3] =  -(2.0/3.0)*s[0][1];
	D[2][3] += 4.0*Ji*W44a*I4a*I4a*a.z*a.z*a.x*a.y;
	D[2][3] += 4.0*Ji*W44b*I4b*I4b*b.z*b.z*b.x*b.y;
	D[2][3] += 4.0*Ji*W44c*I4c*I4c*c.z*c.z*c.x*c.y;
	D[2][3] += 4.0*Ji*W2*(B[2][2]*B[0][1] - B[0][2]*B[1][2]);
	D[2][3] -= (4.0/3.0)*Ji*(W2*(I1*B[0][1] - B2[0][1]) + W44a*I4a*I4a*a.x*a.y + W44b*I4b*I4b*b.x*b.y + W44c*I4c*I4c*c.x*c.y);

	// D[2][4] = 0.5*(c(2,2,1,2) + c(2,2,2,1))
	D[2][4] =  -(2.0/3.0)*s[1][2];
	D[2][4] += 4.0*Ji*W44a*I4a*I4a*a.z*a.z*a.y*a.z;
	D[2][4] += 4.0*Ji*W44b*I4b*I4b*b.z*b.z*b.y*b.z;
	D[2][4] += 4.0*Ji*W44c*I4c*I4c*c.z*c.z*c.y*c.z;
	D[2][4] -= (4.0/3.0)*Ji*(W2*(I1*B[1][2] - B2[1][2]) + W44a*I4a*I4a*a.y*a.z + W44b*I4b*I4b*b.y*b.z + W44c*I4c*I4c*c.y*c.z);

	// D[2][5] = 0.5*(c(2,2,0,2) + c(2,2,2,0))
	D[2][5] =  -(2.0/3.0)*s[0][2];
	D[2][5] += 4.0*Ji*W44a*I4a*I4a*a.z*a.z*a.x*a.z;
	D[2][5] += 4.0*Ji*W44b*I4b*I4b*b.z*b.z*b.x*b.z;
	D[2][5] += 4.0*Ji*W44c*I4c*I4c*c.z*c.z*c.x*c.z;
	D[2][5] -= (4.0/3.0)*Ji*(W2*(I1*B[0][2] - B2[0][2]) + W44a*I4a*I4a*a.x*a.z + W44b*I4b*I4b*b.x*b.z + W44c*I4c*I4c*c.x*c.z);



	// D[3][4] = 0.5*(c(0,1,1,2) + c(0,1,2,1))
	D[3][4] = 4.0*Ji*W44a*I4a*I4a*a.x*a.y*a.y*a.z;
	D[3][4] = 4.0*Ji*W44b*I4b*I4b*b.x*b.y*b.y*b.z;
	D[3][4] = 4.0*Ji*W44c*I4c*I4c*c.x*c.y*c.y*c.z;
	D[3][4] += 2.0*Ji*W2*(B[0][1]*B[1][2] - B[0][2]*B[1][1]);

	// D[3][5] = 0.5*(c(0,1,0,2) + c(0,1,2,0))
	D[3][5] = 4.0*Ji*W44a*I4a*I4a*a.x*a.y*a.x*a.z;
	D[3][5] = 4.0*Ji*W44b*I4b*I4b*b.x*b.y*b.x*b.z;
	D[3][5] = 4.0*Ji*W44c*I4c*I4c*c.x*c.y*c.x*c.z;
	D[3][5] += 2.0*Ji*W2*(B[0][1]*B[0][2] - B[0][0]*B[1][2]);

	// D[4][5] = 0.5*(c(1,2,0,2) + c(1,2,2,0))
	D[4][5] = 4.0*Ji*W44a*I4a*I4a*a.y*a.z*a.x*a.z;
	D[4][5] = 4.0*Ji*W44b*I4b*I4b*b.y*b.z*b.x*b.z;
	D[4][5] = 4.0*Ji*W44c*I4c*I4c*c.y*c.z*c.x*c.z;
	D[4][5] += 2.0*Ji*W2*(B[1][2]*B[0][2] - B[0][1]*B[2][2]);



	// set symmetric components
	D[1][0] = D[0][1]; D[2][0] = D[0][2]; D[3][0] = D[0][3]; D[4][0] = D[0][4]; D[5][0] = D[0][5];
	D[2][1] = D[1][2]; D[3][1] = D[1][3]; D[4][1] = D[1][4]; D[5][1] = D[1][5];
	D[3][2] = D[2][3]; D[4][2] = D[2][4]; D[5][2] = D[2][5];
	D[4][3] = D[3][4]; D[5][3] = D[3][5];
	D[5][4] = D[4][5];

	return tens4ds(D);
}

//-----------------------------------------------------------------------------
// Calculate the deviatoric strain energy density
double FETCNonlinearOrthotropic::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);
    
	// invariants of B
	double I1, I2;
    
	// current local material axis
	vec3d a0, b0, c0, a, b, c;
	double la, lb, lc, lat, lbt, lct;

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber directions
	a0.x = Q[0][0]; b0.x = Q[0][1]; c0.x = Q[0][2];
	a0.y = Q[1][0]; b0.y = Q[1][1]; c0.y = Q[1][2];
	a0.z = Q[2][0]; b0.z = Q[2][1]; c0.z = Q[2][2];
    
	// calculate the current material axes lam*a = F*a0;
	a = F*a0;
	b = F*b0;
	c = F*c0;
    
	// normalize material axis and store fiber stretch
	la  = a.unit();
	lat = la*Jm13; // i.e. lambda tilde
    
	lb  = b.unit();
	lbt = lb*Jm13;
    
	lc  = c.unit();
	lct = lc*Jm13;
    
	// get deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();
    
	// square of B
	mat3ds B2 = B.sqr();
    
	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	I1 = B.tr();
	I2 = 0.5*(I1*I1 - B2.tr());
    
    double sed = m_c1*(I1-3) + m_c2*(I2-3);
    
	// fiber a
	if (lat > 1)
		sed += m_ksi[0]*pow((lat - 1.0), m_beta[0]);
    
	// fiber b
	if (lbt > 1)
		sed += m_ksi[1]*pow((lbt - 1.0), m_beta[1]);
    
	// fiber c
	if (lct > 1)
		sed += m_ksi[2]*pow((lct - 1.0), m_beta[2]);
    
	return sed;
}
