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
#include "FEIsotropicLeeSacks.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEIsotropicLeeSacks, FEElasticMaterial)
	ADD_PARAMETER(m_c0, "c0");
	ADD_PARAMETER(m_c1, "c1");
	ADD_PARAMETER(m_c2, "c2");
	ADD_PARAMETER(m_k , "k");
	ADD_PARAMETER(m_tangent_scale , "tangent_scale");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEIsotropicLeeSacks::FEIsotropicLeeSacks(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_c0 = 0.0;
	m_c1 = 0.0;
	m_c2 = 0.0;
	m_k = 0.0;
	m_tangent_scale = 1.0;
}

//-----------------------------------------------------------------------------
//! Calculates the strain energy density
double FEIsotropicLeeSacks::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate right Cauchy-Green tensor and Lagrange strain tensor
	mat3ds b = pt.LeftCauchyGreen();
	double J = pt.m_J;
	double I1 = b.tr();

	// evaluate volumetric strain energy
	double lnJ = log(J);
	double U = 0.5 * m_k * pow(lnJ, 2);

	// Evaluate exp(Q)
	double Q = m_c2 * (I1 - 3.0) * (I1 - 3.0);
	double eQ = exp(Q);

	// Evaluate the strain energy density
	double sed = m_c0*(0.5*(I1 - 3) - lnJ) + 0.5*m_c1*(eQ - 1.0) + U;

	return sed;
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric stress
mat3ds FEIsotropicLeeSacks::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate deviatoric
	mat3ds b = pt.LeftCauchyGreen();
	double J = pt.m_J;
	double I1 = b.tr();

	// evaluate volumetric strain energy
	double p = m_k * log(J)/J;
	mat3dd I(1.0);

	// Evaluate exp(Q)
	double Q = m_c2 * (I1 - 3.0) * (I1 - 3.0);
	double eQ = exp(Q);

	mat3ds s = I * p + (b - I)*(m_c0/J) + b*(2.0*m_c1*m_c2*(I1 - 3)*eQ/J);

	return s;
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric tangent
tens4ds FEIsotropicLeeSacks::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// apply tangent scale
	double c0 = m_c0 * m_tangent_scale;

	// calculate deviatoric
	mat3ds b = pt.LeftCauchyGreen();
	double J = pt.m_J;
	double I1 = b.tr();

	// we'll need these tensors
	mat3dd I(1.0);
	tens4ds IoI = dyad4s(I);
	tens4ds IxI = dyad1s(I);
	tens4ds bxb = dyad1s(b);

	// evaluate volumetric tangent
	double lnJ = log(J);
	double p = m_k * lnJ / J;
	double pJ = m_k * (1.0 - lnJ) / (J*J);

	tens4ds cp = IxI * (p + J * pJ) + IoI * (2.0 * c0/J -2.0 * p);

	// evaluate strain tangent
	double Q = m_c2 * (I1 - 3.0) * (I1 - 3.0);
	double eQ = exp(Q);

	tens4ds cw = bxb*(4.0*m_c1*m_c2*eQ*(1.0+2.0*m_c2*(I1 - 3.0)*(I1 - 3.0)) / J);

	// put it together
	tens4ds c = cw + cp;

	return c;
}
