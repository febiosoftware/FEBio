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
#include <FECore/ad.h>
#include "adcm.h"

inline double lambdaFromEV(double E, double v)
{
	return v* E / ((1.0 + v) * (1.0 - 2.0 * v));
}

inline double muFromEV(double E, double v)
{
	return 0.5 * E / (1.0 + v);
}

// define the material parameters
BEGIN_FECORE_CLASS(FENeoHookeanAD, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E")->setUnits(UNIT_PRESSURE)->setLongName("Young's modulus");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1, 0.5), "v")->setLongName("Poisson's ratio");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENeoHookeanAD::FENeoHookeanAD(FEModel* pfem) : FEElasticMaterial(pfem) {}

ad::number FENeoHookeanAD::StrainEnergy_AD(FEMaterialPoint& mp, ad::mat3ds& C)
{
	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);
	double lam = lambdaFromEV(E, v);
	double mu = muFromEV(E, v);

	ad::number I1 = C.tr();
	ad::number J = ad::sqrt(C.det());
	ad::number lnJ = ad::log(J);

	ad::number sed = mu * ((I1 - 3) / 2.0 - lnJ) + lam * lnJ * lnJ / 2.0;
	return sed;
}

ad2::number FENeoHookeanAD::StrainEnergy_AD2(FEMaterialPoint& mp, ad2::mat3ds& C)
{
	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);
	double lam = lambdaFromEV(E, v);
	double mu = muFromEV(E, v);

	ad2::number I1 = C.tr();
	ad2::number J = ad2::sqrt(C.det());
	ad2::number lnJ = ad2::log(J);

	ad2::number sed = mu * ((I1 - 3) / 2.0 - lnJ) + lam * lnJ * lnJ / 2.0;
	return sed;
}

ad::mat3ds FENeoHookeanAD::PK2Stress_AD(FEMaterialPoint& mp, ad::mat3ds& C)
{
	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);
	double lam = lambdaFromEV(E, v);
	double mu = muFromEV(E, v);

	ad::mat3ds I(1.0);
	ad::mat3ds Ci = C.inverse();
	ad::number J = ad::sqrt(C.det());
	ad::number lnJ = ad::log(J);

	ad::mat3ds S = (I - Ci) * mu + Ci * (lam * lnJ);

	return S;
}

mat3ds FENeoHookeanAD::Stress(FEMaterialPoint& mp)
{
	// calculate PK2 stress
//	mat3ds S = ad::PK2Stress<FENeoHookeanAD>(this, mp);
	mat3ds S = ad2::PK2Stress<FENeoHookeanAD>(this, mp);

	// push-forward to obtain Cauchy-stress
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3ds s = pt.push_forward(S);

	return s;
}

tens4ds FENeoHookeanAD::Tangent(FEMaterialPoint& mp)
{
	// calculate material tangent
//	tens4ds C4 = ad::Tangent<FENeoHookeanAD>(this, mp);
	tens4ds C4 = ad2::Tangent<FENeoHookeanAD>(this, mp);

	// push forward to get spatial tangent
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	tens4ds c4 = pt.push_forward(C4);

	return c4;
}

double FENeoHookeanAD::StrainEnergyDensity(FEMaterialPoint& mp)
{
	return ad::StrainEnergy<FENeoHookeanAD>(this, mp);
}

mat3ds FENeoHookeanAD::PK2Stress(FEMaterialPoint& mp, const mat3ds ES)
{
	mat3ds C = mat3dd(1) + ES * 2;
	mat3ds S = ad::PK2Stress<FENeoHookeanAD>(this, mp, C);
	return S;
}

tens4dmm FENeoHookeanAD::MaterialTangent(FEMaterialPoint& mp, const mat3ds ES)
{
	mat3ds C = mat3dd(1) + ES * 2;
	tens4ds C4 = ad::Tangent<FENeoHookeanAD>(this, mp, C);
	return tens4dmm(C4);
}
