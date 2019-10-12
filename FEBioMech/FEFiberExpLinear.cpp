/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEFiberExpLinear.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEFiberExpLinear, FEElasticFiberMaterial)
	ADD_PARAMETER(m_c3  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c3");
	ADD_PARAMETER(m_c4  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c4");
	ADD_PARAMETER(m_c5  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c5");
	ADD_PARAMETER(m_lam1, FE_RANGE_GREATER_OR_EQUAL(1.0), "lambda");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFiberExpLinear::FEFiberExpLinear(FEModel* pfem) : FEElasticFiberMaterial(pfem)
{
	m_c3 = 0;
	m_c4 = 0;
	m_c5 = 0;
	m_lam1 = 1.0;
}

//-----------------------------------------------------------------------------
//! Calculate the fiber stress
mat3ds FEFiberExpLinear::FiberStress(FEMaterialPoint& mp, const vec3d& a0)
{
	// get the material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the spatial fiber vector and stretch
	vec3d a = pt.m_F*a0;
	double l = a.unit();

	// other stuff we need
	mat3ds A = dyad(a);
	double J = pt.m_J;

	// fiber stress
	mat3ds s; s.zero();

	// calculate fiber stress
	if (l >= 1.0)
	{
		double Wl = 0.0;
		if (l < m_lam1)
		{
			Wl = m_c3*(exp(m_c4*(l - 1.0)) - 1.0);
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1 - 1.0)) - 1.0) - m_c5*m_lam1;
			Wl = m_c5*l + c6;
		}
		s += A*(Wl / J);
	}

	return s;
}

//-----------------------------------------------------------------------------
//! Calculate the fiber tangent
tens4ds FEFiberExpLinear::FiberTangent(FEMaterialPoint& mp, const vec3d& a0)
{
	// get material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the spatial fiber axis
	vec3d a = pt.m_F*a0;
	double l = a.unit();

	// Invariants of B (= invariants of C)
	double J = pt.m_J;
	double I4 = l*l;

	// some useful tensors
	mat3dd I(1.0);
	mat3ds A = dyad(a);
	tens4ds IxI = dyad1s(I);
	tens4ds IoI = dyad4s(I);
	tens4ds AxA = dyad1s(A);

	// fiber tangent
	tens4ds c(0.0);
	if (l >= 1.0)
	{
		double Fl = 0.0, Fll = 0.0;
		if (l < m_lam1)
		{
			Fl = m_c3*(exp(m_c4*(l - 1.0)) - 1.0) / l;
			Fll = -m_c3*(exp(m_c4*(l - 1.0)) - 1.0) / (l*l) + m_c3*m_c4*exp(m_c4*(l - 1.0)) / l;
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1 - 1.0)) - 1.0) - m_c5*m_lam1;
			Fl = m_c5 + c6 / l;
			Fll = -c6 / (l*l);
		}

		double W44 = (Fll - Fl / l) / (4 * l*l);

		c += AxA*(4.0*W44*I4*I4 / J);
	}

	return c;
}

//-----------------------------------------------------------------------------
//! Calculate the fiber strain energy density
double FEFiberExpLinear::FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0)
{
	// TODO: implement this
	return 0;
}
