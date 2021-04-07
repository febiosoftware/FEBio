/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
BEGIN_FECORE_CLASS(FEUncoupledFiberExpLinear, FEElasticFiberMaterialUC);
	ADD_PARAMETER(m_c3  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c3");
	ADD_PARAMETER(m_c4  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c4");
	ADD_PARAMETER(m_c5  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c5");
	ADD_PARAMETER(m_lam1, FE_RANGE_GREATER_OR_EQUAL(1.0), "lambda");
	ADD_PARAMETER(m_fiber, "fiber");
	ADD_PARAMETER(m_epsf, "epsilon_scale");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEUncoupledFiberExpLinear::FEUncoupledFiberExpLinear(FEModel* pfem) : FEElasticFiberMaterialUC(pfem)
{
	m_c3 = m_c4 = m_c5 = 0;
	m_lam1 = 1;

	m_epsf = 0.0;
}

//-----------------------------------------------------------------------------
//! Fiber material stress
mat3ds FEUncoupledFiberExpLinear::DevFiberStress(FEMaterialPoint &mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Ji = 1.0 / J;
	double Jm13 = pow(J, -1.0 / 3.0);
	double twoJi = 2.0*Ji;

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// invariant I4
	double I4 = lamd*lamd;

	// strain energy derivative
	double W4 = 0;
	if (lamd >= 1)
	{
		double lamdi = 1.0 / lamd;
		double Wl;
		if (lamd < m_lam1)
		{
			Wl = lamdi*m_c3*(exp(m_c4*(lamd - 1)) - 1);
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1 - 1)) - 1) - m_c5*m_lam1;
			Wl = lamdi*(m_c5*lamd + c6);
		}
		W4 = 0.5*lamdi*Wl;
	}
	else
	{
		W4 = 0;
	}

	// calculate dyad of a: AxA = (a x a)
	mat3ds AxA = dyad(a);

	// ---
	// calculate FdWf/dCFt = I4*W4*(a x a)
	mat3ds T = AxA*(W4*I4);

	// calculate stress: 
	mat3ds s = T.dev()*twoJi;

	return s;
}

//-----------------------------------------------------------------------------
//! Fiber material tangent
tens4ds FEUncoupledFiberExpLinear::DevFiberTangent(FEMaterialPoint &mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);
	double Ji = 1.0 / J;

	// calculate current local material axis
	vec3d a = F*a0;

	double lam = a.unit();

	// deviatoric stretch
	double lamd = lam*Jm13;

	double I4 = lamd*lamd;

	const double eps = m_epsf*std::numeric_limits<double>::epsilon();

	double W4, W44;
	if (lamd >= 1 + eps)
	{
		double lamdi = 1.0 / lamd;
		double Wl, Wll;
		if (lamd < m_lam1)
		{
			Wl = lamdi*m_c3*(exp(m_c4*(lamd - 1)) - 1);
			Wll = m_c3*lamdi*(m_c4*exp(m_c4*(lamd - 1)) - lamdi*(exp(m_c4*(lamd - 1)) - 1));
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1 - 1)) - 1) - m_c5*m_lam1;
			Wl = lamdi*(m_c5*lamd + c6);
			Wll = -c6*lamdi*lamdi;
		}
		W4 = 0.5*lamdi*Wl;
		W44 = 0.25*lamdi*lamdi*(Wll - lamdi*Wl);
	}
	else
	{
		W4 = 0;
		W44 = 0;
	}

	// --- calculate tangent ---

	// calculate dWdC:C
	double WC = W4*I4;

	// calculate C:d2WdCdC:C
	double CWWC = W44*I4*I4;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds Id4 = dyad4s(I);

	mat3ds AxA = dyad(a);
	tens4ds AxAxAxA = dyad1s(AxA);

	tens4ds cw = AxAxAxA*(4.0*Ji*W44*I4*I4) - dyad1s(I, AxA)*(4.0 / 3.0*Ji*W44*I4*I4);

	tens4ds c = (Id4 - IxI / 3.0)*(4.0 / 3.0*Ji*WC) + IxI*(4.0 / 9.0*Ji*CWWC) + cw;

	return c;
}

//-----------------------------------------------------------------------------
//! Fiber material strain energy density
double FEUncoupledFiberExpLinear::DevFiberStrainEnergyDensity(FEMaterialPoint &mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// strain energy density
	double sed = 0.0;
	if (lamd >= 1)
	{
		if (lamd < m_lam1)
		{
			sed = m_c3*(exp(-m_c4)*
				(expint_Ei(m_c4*lamd) - expint_Ei(m_c4))
				- log(lamd));
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1 - 1)) - 1) - m_c5*m_lam1;
			sed = m_c5*(lamd - 1) + c6*log(lamd);
		}
	}
	// --- active contraction contribution to sed is zero ---

	return sed;
}
