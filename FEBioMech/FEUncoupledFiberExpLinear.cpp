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
#include "FEUncoupledFiberExpLinear.h"
#include <stdlib.h>
#include <limits>
#include <FECore/expint_Ei.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFiberExpLinearUC, FEFiberMaterialUncoupled);
ADD_PARAMETER(m_c3, FE_RANGE_GREATER_OR_EQUAL(0.0), "c3");
ADD_PARAMETER(m_c4, FE_RANGE_GREATER_OR_EQUAL(0.0), "c4");
ADD_PARAMETER(m_c5, FE_RANGE_GREATER_OR_EQUAL(0.0), "c5");
ADD_PARAMETER(m_lam1, FE_RANGE_GREATER_OR_EQUAL(1.0), "lambda");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFiberExpLinearUC::FEFiberExpLinearUC(FEModel* pfem) : FEFiberMaterialUncoupled(pfem)
{
	m_c3 = 0;
	m_c4 = 0;
	m_c5 = 0;
	m_lam1 = 1;
}

//-----------------------------------------------------------------------------
//! Fiber material stress
mat3ds FEFiberExpLinearUC::DevFiberStress(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Ji = 1.0 / J;
	double Jm13 = pow(J, -1.0 / 3.0);
	double twoJi = 2.0 * Ji;

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F * a0;

	// normalize material axis and store fiber stretch
	double lam = a.unit();
	double lamd = lam * Jm13; // i.e. lambda tilde

	// invariant I4
	double I4 = lamd * lamd;

	// calculate stress:
	mat3ds s; s.zero();

	if (lamd >= 1)
	{
		double c3 = m_c3(mp);
		double c4 = m_c4(mp);
		double c5 = m_c5(mp);
		double lam1 = m_lam1(mp);

		// calculate dyad of a: AxA = (a x a)
		mat3ds AxA = dyad(a);

		if (c3 == 0) {
			c3 = c5 / c4 * exp(-c4 * (lam1 - 1));
		}

		// calculate fiber stress
		double sn = 0.0;
		if (lamd < lam1)
		{
			sn = c3 * (exp(c4 * (lamd - 1.0)) - 1.0);
		}
		else
		{
			double c6 = c3 * (exp(c4 * (lam1 - 1)) - 1) - c5 * lam1;
			sn = c5 * lamd + c6;
		}
		mat3ds T = AxA * (sn / J);
		s = T.dev();
	}

	return s;
}

//-----------------------------------------------------------------------------
//! Fiber material tangent
tens4ds FEFiberExpLinearUC::DevFiberTangent(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);
	double Ji = 1.0 / J;

	// calculate current local material axis
	vec3d a = F * a0;

	double lam = a.unit();

	// deviatoric stretch
	double lamd = lam * Jm13;

	double I4 = lamd * lamd;

	const double eps = 0;// std::numeric_limits<double>::epsilon();

	tens4ds c; c.zero();

	if (lamd >= 1 + eps)
	{
		double c3 = m_c3(mp);
		double c4 = m_c4(mp);
		double c5 = m_c5(mp);
		double lam1 = m_lam1(mp);

		mat3dd I(1);    // Identity
		tens4ds IxI = dyad1s(I);
		tens4ds Id4 = dyad4s(I);

		mat3ds AxA = dyad(a);
		tens4ds AxAxAxA = dyad1s(AxA);

		if (c3 == 0) {
			c3 = c5 / c4 * exp(-c4 * (lam1 - 1));
		}

		double sn = 0;
		double cn = 0;
		if (lamd < lam1) {
			sn = c3 * (exp(c4 * (lamd - 1.0)) - 1.0);
			cn = c3 * (2 + exp(c4 * (lamd - 1)) * (c4 * lamd - 2));
		}
		else {
			double c6 = c3 * (exp(c4 * (lam1 - 1)) - 1) - c5 * lam1;
			sn = c5 * lamd + c6;
			cn = -c5 * lamd - 2 * c6;
		}
		mat3ds T = AxA * (sn / J);
		c = AxAxAxA * (cn / J);
		c += -1. / 3. * (ddots(c, IxI) - IxI * (c.tr() / 3.))
			+ 2. / 3. * ((Id4 - IxI / 3.) * T.tr() - dyad1s(T.dev(), I));
	}

	return c;
}

//-----------------------------------------------------------------------------
//! Fiber material strain energy density
double FEFiberExpLinearUC::DevFiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F * a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam * Jm13; // i.e. lambda tilde

	// strain energy density
	double sed = 0.0;
	if (lamd >= 1)
	{
		double c3 = m_c3(mp);
		double c4 = m_c4(mp);
		double c5 = m_c5(mp);
		double lam1 = m_lam1(mp);

		if (c3 == 0) c3 = c5 / c4 * exp(-c4 * (lam1 - 1));

		if (lamd < lam1)
		{
			sed = c3 * exp(-c4) * (expint_Ei(c4 * lamd) - expint_Ei(c4)) - c3 * log(lamd);
		}
		else
		{
			double c6 = c3 * (exp(c4 * (lam1 - 1)) - 1) - c5 * lam1;
			sed = c5 * (lamd - lam1) + c6 * log(lamd / lam1)
				+ c3 * exp(-c4) * (expint_Ei(c4 * lam1) - expint_Ei(c4)) - c3 * log(lam1);
		}
	}

	return sed;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEUncoupledFiberExpLinear, FEElasticFiberMaterialUC);
ADD_PARAMETER(m_fib.m_c3, FE_RANGE_GREATER_OR_EQUAL(0.0), "c3");
ADD_PARAMETER(m_fib.m_c4, FE_RANGE_GREATER_OR_EQUAL(0.0), "c4");
ADD_PARAMETER(m_fib.m_c5, FE_RANGE_GREATER_OR_EQUAL(0.0), "c5");
ADD_PARAMETER(m_fib.m_lam1, FE_RANGE_GREATER_OR_EQUAL(1.0), "lambda");
END_FECORE_CLASS();
