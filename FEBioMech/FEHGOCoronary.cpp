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
#include "FEHGOCoronary.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEHGOCoronary, FEUncoupledMaterial)
	ADD_PARAMETER(m_rho, "rho");
	ADD_PARAMETER(m_k1 , "k1");
	ADD_PARAMETER(m_k2 , "k2");

	ADD_PARAMETER(m_fiber, "fiber");
	
END_FECORE_CLASS();

//////////////////////////////////////////////////////////////////////
// FETransIsoMooneyRivlin
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
FEHGOCoronary::FEHGOCoronary(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	m_rho = 0.0;
	m_k1 = 0.0;
	m_k2 = 1.0;
	m_fiber = vec3d(1,0,0);
}

//-----------------------------------------------------------------------------
mat3ds FEHGOCoronary::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// fiber vector
	mat3d Q = GetLocalCS(mp);
	vec3d a0 = Q*m_fiber(mp); a0.unit();
	vec3d a = F * a0;
	double lam = Jm13 * a.unit();
	mat3ds m = dyad(a);

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I4 = lam * lam;

	// material parameters
	double rho = m_rho;
	double k1 = m_k1;
	double k2 = m_k2;

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double w = k2 * ((1.0 - rho) * (I1 - 3.0) * (I1 - 3.0) + rho*(I4 - 1.0) * (I4 - 1.0));
	double w1  = 2.0 * k2 * (1.0 - rho) * (I1 - 3.0);
	double w4  = 2.0 * k2 * rho * (I4 - 1.0);
	double kew = k1*exp(w) / k2;

	double W1 = kew * w1;
	double W4 = kew * w4;
	// ------------------------------------------------

	// calculate T = F*dW/dC*Ft
	mat3ds T = B * W1 + m*(W4*I4);

	// calculate stress s = pI + 2/J * dev(T) 
	mat3ds s = T.dev() * (2.0 / J);

	return s;
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4ds FEHGOCoronary::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Ji = 1.0 / pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// fiber vector
	mat3d Q = GetLocalCS(mp);
	vec3d a0 = Q*m_fiber(mp); a0.unit();
	vec3d a = F * a0;
	double lam = Jm13 * a.unit();
	mat3ds m = dyad(a);

	// Invariants of B (= invariants of C)
	double J1 = B.tr();
	double J4 = lam * lam;

	// material parameters
	double rho = m_rho;
	double k1 = m_k1;
	double k2 = m_k2;

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double w = k2 * ((1.0 - rho) * (J1 - 3.0) * (J1 - 3.0) + rho * (J4 - 1.0) * (J4 - 1.0));
	double w1 = 2.0 * k2 * (1.0 - rho) * (J1 - 3.0);
	double w11 = 2.0 * k2 * (1.0 - rho);
	double w4 = 2.0 * k2 * rho * (J4 - 1.0);
	double w44 = 2.0 * k2 * rho;
	double kew = k1 * exp(w) / k2;

	double W1 = kew * w1;
	double W4 = kew * w4;
	double W11 = kew * (w1*w1 + w11);
	double W44 = kew * (w4*w4 + w44);
	double W14 = kew * (w1*w4);
	// ------------------------------------

	// calculate T = F*dW/dC*Ft
	mat3ds T = B * W1 + m * (W4 * J4);

	// calculate stress s = pI + 2/J * dev(T) 
	mat3ds devs = T.dev() * (2.0 / J);

	// calculate dWdC:C
	double WC = W1 * J1 + W4 * J4;

	// calculate C:d2WdCdC:C
	double CWWC = W11 * J1 * J1 + 2.0 * W14 * J1 * J4 + W44 * J4 * J4;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4 = dyad4s(B);
	tens4ds Bxm = dyad1s(B, m);
	tens4ds mxm = dyad1s(m);

	// d2W/dCdC:C
	mat3ds WCCxC = B * (W11 * J1) + B*(W14*J4) + m*(W14*J4*J1) + m*(W44*J4*J4);

	tens4ds Cw = BxB * W11 + Bxm * (W14 * J4) + mxm * (W44 * J4 * J4);

	tens4ds cw = Cw*(4.0 *Ji) - dyad1s(WCCxC, I) * (4.0 / 3.0 * Ji) + IxI * (4.0 / 9.0 * Ji * CWWC);
	tens4ds c = dyad1s(devs, I) * (-2.0 / 3.0) + (I4 - IxI / 3.0) * (4.0 / 3.0 * Ji * WC) + cw;

	return c;
}

//-----------------------------------------------------------------------------
double FEHGOCoronary::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// fiber vector
	mat3d Q = GetLocalCS(mp);
	vec3d a0 = Q*m_fiber(mp); a0.unit();
	vec3d a = F * a0;
	double lam = Jm13 * a.unit();

	// Invariants of B (= invariants of C)
	double J1 = B.tr();
	double J4 = lam * lam;

	// material parameters
	double rho = m_rho;
	double k1 = m_k1;
	double k2 = m_k2;

	// calculate sed
	double w = k2 * ((1.0 - rho) * (J1 - 3.0) * (J1 - 3.0) + rho * (J4 - 1.0) * (J4 - 1.0));
	double sed = (k1/k2)*(exp(w) - 1.0);

	return sed;
}
