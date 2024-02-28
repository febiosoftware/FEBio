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
#include "FENeoHookeanAD.h"
#include <FECore/fad.h>

double _lam = 0.0;
double _mu = 0.0;

// the strain energy function to use with AD
fad::number W_nh(const fad::Mat3ds& C)
{
	fad::number I1 = C.Trace();
	fad::number J2 = C.Det();
	fad::number J = fad::Sqrt(J2);
	fad::number lnJ = fad::Ln(J);

	fad::number sed = _mu * ((I1 - 3) / 2.0 - lnJ) + _lam * lnJ * lnJ / 2.0;
	return sed;
}

// define the material parameters
BEGIN_FECORE_CLASS(FENeoHookeanAD, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E")->setUnits(UNIT_PRESSURE)->setLongName("Young's modulus");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1, 0.5), "v")->setLongName("Poisson's ratio");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENeoHookeanAD::FENeoHookeanAD(FEModel* pfem) : FEElasticMaterial(pfem) {}

//-----------------------------------------------------------------------------
mat3ds FENeoHookeanAD::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);

	// lame parameters
	_lam = v * E / ((1 + v) * (1 - 2 * v));
	_mu = 0.5 * E / (1 + v);

	// calculate PK2 stress
	mat3ds C = pt.RightCauchyGreen();
	mat3ds S = fad::Derive(W_nh, C)*2.0;

	// push-forward to obtain Cauchy-stress
	mat3ds s = pt.push_forward(S);

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FENeoHookeanAD::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	double detF = pt.m_J;

	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);

	// lame parameters
	double lam = v * E / ((1 + v) * (1 - 2 * v));
	double mu = 0.5 * E / (1 + v);

	double lam1 = lam / detF;
	double mu1 = (mu - lam * log(detF)) / detF;

	mat3dd I(1);

	return dyad1s(I) * lam1 + dyad4s(I) * (2 * mu1);
}

//-----------------------------------------------------------------------------
double FENeoHookeanAD::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);

	// lame parameters
	_lam = v * E / ((1 + v) * (1 - 2 * v));
	_mu = 0.5 * E / (1 + v);

	mat3ds C = pt.RightCauchyGreen();
	double sed = fad::Evaluate(W_nh, C);

	return sed;
}

//-----------------------------------------------------------------------------
mat3ds FENeoHookeanAD::PK2Stress(FEMaterialPoint& pt, const mat3ds ES)
{
	// Identity
	mat3dd I(1);

	// calculate right Cauchy-Green tensor
	mat3ds C = I + ES * 2;
	mat3ds Ci = C.inverse();

	double detF = sqrt(C.det());
	double lndetF = log(detF);

	// get the material parameters
	double E = m_E(pt);
	double v = m_v(pt);

	// lame parameters
	double lam = v * E / ((1 + v) * (1 - 2 * v));
	double mu = 0.5 * E / (1 + v);

	// calculate stress
	mat3ds S = (I - Ci) * mu + Ci * (lam * lndetF);

	return S;
}

//-----------------------------------------------------------------------------
tens4dmm FENeoHookeanAD::MaterialTangent(FEMaterialPoint& pt, const mat3ds ES)
{
	// calculate right Cauchy-Green tensor
	mat3ds C = mat3dd(1) + ES * 2;
	mat3ds Ci = C.inverse();
	double J = sqrt(C.det());

	// get the material parameters
	double E = m_E(pt);
	double v = m_v(pt);

	// lame parameters
	double lam = v * E / ((1 + v) * (1 - 2 * v));
	double mu = 0.5 * E / (1 + v);

	tens4dmm c = dyad1s(Ci) * lam + dyad4s(Ci) * (2 * (mu - lam * log(J)));

	return c;
}
