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
#include "FEPolynomialHyperElastic.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPolynomialHyperElastic, FEElasticMaterial)
//	ADD_PARAMETER(m_c[0][0], "c00");  // should always be zero
	ADD_PARAMETER(m_c[0][1], "c01");
	ADD_PARAMETER(m_c[0][2], "c02");
	ADD_PARAMETER(m_c[1][0], "c10");
	ADD_PARAMETER(m_c[1][1], "c11");
	ADD_PARAMETER(m_c[1][2], "c12");
	ADD_PARAMETER(m_c[2][0], "c20");
	ADD_PARAMETER(m_c[2][1], "c21");
	ADD_PARAMETER(m_c[2][2], "c22");
	ADD_PARAMETER(m_D1, "D1");
	ADD_PARAMETER(m_D2, "D2");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPolynomialHyperElastic::FEPolynomialHyperElastic(FEModel* pfem) : FEElasticMaterial(pfem) 
{
	m_D1 = m_D2 = 0.0;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j) m_c[i][j] = 0.0;
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric stress
mat3ds FEPolynomialHyperElastic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5 * (I1 * I1 - B2.tr());

	// get material parameters
	double c[3][3] = { 0 };
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			c[i][j] = m_c[i][j](mp);
	c[0][0] = 0.0;

	// --- put strain energy derivatives here ---
	// Wi = dW/dIi
	double I1p1 = I1 - 3.0;
	double I1p2 = I1p1 * I1p1;
	double I2p1 = I2 - 3.0;
	double I2p2 = I2p1 * I2p1;

	double W1 = c[1][0] + c[1][1] * I2p1 + c[1][2] * I2p2 + 2.0 * I1p1 * (c[2][0] + c[2][1] * I2p1 + c[2][2] * I2p2);
	double W2 = c[0][1] + c[1][1] * I1p1 + c[2][1] * I1p2 + 2.0 * I2p1 * (c[0][2] + c[1][2] * I1p1 + c[2][2] * I1p2);
	// ---

	// calculate T = F*dW/dC*Ft
	// T = F*dW/dC*Ft
	mat3ds T = B * (W1 + W2 * I1) - B2 * W2;

	mat3ds devs = T.dev() * (2.0 / J);

	pt.m_p = UJ(pt.m_J);

	return mat3dd(pt.m_p) + devs;
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric tangent
tens4ds FEPolynomialHyperElastic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;
	double Ji = 1.0 / J;

	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double I2 = 0.5 * (I1 * I1 - B2.tr());

	// get material parameters
	double c[3][3] = { 0 };
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			c[i][j] = m_c[i][j](mp);
	c[0][0] = 0.0;

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double I1p1 = I1 - 3.0;
	double I1p2 = I1p1 * I1p1;
	double I2p1 = I2 - 3.0;
	double I2p2 = I2p1 * I2p1;

	double W1 = c[1][0] + c[1][1] * I2p1 + c[1][2] * I2p2 + 2.0 * I1p1 * (c[2][0] + c[2][1] * I2p1 + c[2][2] * I2p2);
	double W2 = c[0][1] + c[1][1] * I1p1 + c[2][1] * I1p2 + 2.0 * I2p1 * (c[0][2] + c[1][2] * I1p1 + c[2][2] * I1p2);

	double W11 = 2.0 * (c[2][0] + c[2][1] * I2p1 + c[2][2] * I2p2);
	double W22 = 2.0 * (c[0][2] + c[1][2] * I1p1 + c[2][2] * I1p2);

	double W12 = c[1][1] + 2.0 * c[1][2] * I2p1 + 2.0 * c[2][1] * I1p1 + 4.0 * c[2][2]*I1p1 * I2p1;
	double W21 = W12;
	// ---

	// define some helper constants
	double g1 = W11 + 2.0 * W12 * I1 + W22 * I1 * I1 + W2;
	double g2 = W12 + I1 * W22;
	double I2b = B2.tr();

	// calculate dWdC:C
	double WC = W1 * I1 + 2 * W2 * I2;

	// deviatoric cauchy-stress
	mat3ds T = B * (W1 + W2 * I1) - B2 * W2;
	mat3ds devs = T.dev() * (2.0 / J);

	// Identity tensor
	mat3ds I(1, 1, 1, 0, 0, 0);

	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B2xB2 = dyad1s(B2);
	tens4ds B4 = dyad4s(B);

	// calculate chat
	tens4ds ch = BxB * g1 - dyad1s(B, B2) * g2 + B2xB2 * W22 - B4 * W2;

	// calculate I:ch:I
	double IchI = W11*I1*I1 + 2.0*W12*I1*I2 + 2.0*W21*I1*I2 + 4.0*W22*I2*I2 + 2.0*W2*I2;

	// 1:ch
	mat3ds Ich = B * (g1*I1-g2*I2b) + B2 * (W22*I2b - g2*I1 - W2);

	tens4ds cw = ch - dyad1s(Ich, I) * (1.0 / 3.0) + IxI * (1.0 / 9.0 * IchI);

	tens4ds cs = dyad1s(devs, I)* (-2.0 / 3.0) + (I4 - IxI / 3.0) * (4.0 / 3.0 * Ji * WC);
	tens4ds ce = cs + cw*(4.0 * Ji);

	return ce + (IxI - I4 * 2) * pt.m_p + IxI * (UJJ(pt.m_J) * pt.m_J);
}

//-----------------------------------------------------------------------------
//! calculate deviatoric strain energy density
double FEPolynomialHyperElastic::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get material parameters
	double c[3][3] = { 0 };
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			c[i][j] = m_c[i][j](mp);

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5 * (I1 * I1 - B2.tr());

	double W = 0.0;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			W += c[i][j] * pow(I1 - 3, i) * pow(I2 - 3, j);

	return U(pt.m_J) + W;
}

double FEPolynomialHyperElastic::U(double J)
{
	return m_D1 * pow(J - 1.0, 2.0) + m_D2 * pow(J - 1.0, 4.0);
}

double FEPolynomialHyperElastic::UJ(double J)
{
	return 2.0*m_D1 * (J - 1.0) + 4.0*m_D2 * pow(J - 1.0, 3.0);
}

double FEPolynomialHyperElastic::UJJ(double J)
{
	return 2.0*m_D1 + 12.0*m_D2 * pow(J - 1.0, 2.0);
}
