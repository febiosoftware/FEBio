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
#include "FEIsotropicLeeSacksUncoupled.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEIsotropicLeeSacksUncoupled, FEUncoupledMaterial)
	ADD_PARAMETER(m_c0, "c0");
	ADD_PARAMETER(m_c1, "c1");
	ADD_PARAMETER(m_c2, "c2");
	ADD_PARAMETER(m_tangent_scale, "tangent_scale");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEIsotropicLeeSacksUncoupled::FEIsotropicLeeSacksUncoupled(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	m_c0 = 0.0;
	m_c1 = 0.0;
	m_c2 = 0.0;
	m_tangent_scale = 1.0;
}

//-----------------------------------------------------------------------------
//! Calculates the strain energy density
double FEIsotropicLeeSacksUncoupled::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate left Cauchy-Green tensor and invariants
	mat3ds b = pt.DevLeftCauchyGreen();
	double I1 = b.tr();

	// Evaluate exp(Q)
	double Q = m_c2 * (I1 - 3.0) * (I1 - 3.0);
	double eQ = exp(Q);

	// Evaluate the strain energy density
	double sed = 0.5*m_c0*(I1 - 3.0) + 0.5*m_c1*(eQ - 1.0);

	return sed;
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric stress
mat3ds FEIsotropicLeeSacksUncoupled::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate left Cauchy-Green tensor and invariants
	mat3ds b = pt.DevLeftCauchyGreen();
	double I1 = b.tr();
	double J = pt.m_J;

	// Evaluate exp(Q)
	double Q = m_c2 * (I1 - 3.0) * (I1 - 3.0);
	double eQ = exp(Q);

	// evaluate deviatoric stress
	mat3ds T = b * (m_c0 + 2.0 * m_c1 * m_c2 * (I1 - 3) * eQ);
	mat3ds s = T.dev() / J;
	return s;
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric tangent
tens4ds FEIsotropicLeeSacksUncoupled::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// apply tangent scale
	double c0 = m_c0 * m_tangent_scale;

	// calculate left Cauchy-Green tensor and invariants
	mat3ds b = pt.DevLeftCauchyGreen();
	double I1 = b.tr();
	double J = pt.m_J;

	// Evaluate exp(Q)
	double Q = m_c2 * (I1 - 3.0) * (I1 - 3.0);
	double eQ = exp(Q);
	double g = 4.0 * m_c1 * m_c2 * eQ * (1.0 + 2.0 * m_c2 * (I1 - 3.0) * (I1 - 3.0));

	// we'll need these tensors
	mat3dd I(1.0);
	tens4ds IoI = dyad4s(I);
	tens4ds IxI = dyad1s(I);
	tens4ds bxb = dyad1s(b);

	// evaluate stress tensor T = F*S*F^t
	mat3ds T = b * (c0 + 2.0 * m_c1 * m_c2 * (I1 - 3) * eQ);

	// evaluate tangents
	tens4ds cT = (IoI - IxI / 3.0) * (2.0 * T.tr() / 3.0) - dyad1s(T.dev(), I)*(2.0/3.0);

	tens4ds cw = bxb * g;
	mat3ds chat = b * (g * I1);
	double cbar = g * I1 * I1;

	tens4ds c = cT + cw - dyad1s(chat, I) / 3.0 + IxI*(cbar / 9.0);

	return c / J;
}
