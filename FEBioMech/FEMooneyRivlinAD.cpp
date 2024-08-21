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
#include "FEMooneyRivlinAD.h"
#include "adcm.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEMooneyRivlinAD, FEUncoupledMaterial)
	ADD_PARAMETER(m_c1, "c1")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c2, "c2")->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Calculate the deviatoric stress
mat3ds FEMooneyRivlinAD::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get material parameters
	double c1 = m_c1(mp);
	double c2 = m_c2(mp);

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

	// --- TODO: put strain energy derivatives here ---
	//
	// W = C1*(I1 - 3) + C2*(I2 - 3)
	//
	// Wi = dW/dIi
	double W1 = c1;
	double W2 = c2;
	// ---

	// calculate T = F*dW/dC*Ft
	// T = F*dW/dC*Ft
	mat3ds T = B * (W1 + W2 * I1) - B2 * W2;

	return T.dev() * (2.0 / J);
}

ad::mat3ds FEMooneyRivlinAD::PK2Stress_AD(FEMaterialPoint& mp, ad::mat3ds& C)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get material parameters
	double c1 = m_c1(mp);
	double c2 = m_c2(mp);

	// determinant of deformation gradient
	double J = pt.m_J;
	double Jm23 = pow(J, -1.0 / 3.0);

	ad::mat3ds C2 = C.sqr();
	ad::mat3ds Ci = C.inverse();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	ad::number I1 = C.tr();
	ad::number I2 = 0.5 * (I1 * I1 - C2.tr());

	// calculate T = dW/dC
	ad::mat3ds I(1.0);
	ad::mat3ds T = I*c1 + C*(0.5*c2);

	// calculte S = 2*DEV[T]
	ad::mat3ds S = (T - Ci * (T.dotdot(C) / 3.0))*(2.0*Jm23);

	return S;
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric tangent
tens4ds FEMooneyRivlinAD::DevTangent(FEMaterialPoint& mp)
{
	// calculate material tangent
	tens4ds C4 = ad::Tangent<FEMooneyRivlinAD>(this, mp);

	// push forward to get spatial tangent
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	tens4ds c4 = pt.push_forward(C4);

	return c4;
}

//-----------------------------------------------------------------------------
//! calculate deviatoric strain energy density
double FEMooneyRivlinAD::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get material parameters
	double c1 = m_c1(mp);
	double c2 = m_c2(mp);

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5 * (I1 * I1 - B2.tr());

	//
	// W = C1*(I1 - 3) + C2*(I2 - 3)
	//
	double sed = c1 * (I1 - 3) + c2 * (I2 - 3);

	return sed;
}
