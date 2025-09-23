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
#include <limits>
#include "FEFiberPow.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberPow, FEFiberMaterial)
	ADD_PARAMETER(m_ksi, FE_RANGE_GREATER(0.0), "ksi")->setUnits(UNIT_PRESSURE)->setLongName("fiber modulus");
	ADD_PARAMETER(m_beta, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta")->setLongName("power exponent");
	ADD_PARAMETER(m_tension_only, "tension_only");
END_FECORE_CLASS();

FEFiberPow::FEFiberPow(FEModel* pfem) : FEFiberMaterial(pfem)
{
	m_ksi = 0.0;
	m_beta = 2.0;
	m_epsf = 0.0;
	m_tension_only = true;
}

mat3ds FEFiberPow::FiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// initialize material constants
	double ksi = m_ksi(mp);
	double beta = m_beta(mp);

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;

	// loop over all integration points
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds s;

	// Calculate In
	double In = n0 * (C * n0);

	// only take fibers in tension into consideration
	if (!m_tension_only || (In - 1 >= eps))
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F * n0 / sqrt(In);

		// calculate the outer product of nt
		mat3ds N = dyad(nt);

		// calculate the fiber stress magnitude
		double sn = 2 * In * ksi * beta * pow(In - 1, beta - 1);

		// calculate the fiber stress
		s = N * (sn / J);
	}
	else
	{
		s.zero();
	}

	return s;
}

tens4ds FEFiberPow::FiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// initialize material constants
	double ksi = m_ksi(mp);
	double beta = m_beta(mp);

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;

	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
	tens4ds c;

	// Calculate In
	double In = n0 * (C * n0);

	// only take fibers in tension into consideration
	const double eps = m_epsf * std::numeric_limits<double>::epsilon();
	if (!m_tension_only || (In >= 1 + eps))
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F * n0 / sqrt(In);

		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		tens4ds NxN = dyad1s(N);

		// calculate modulus
		double cn = 4 * In * In * ksi * beta * (beta - 1) * pow(In - 1, beta - 2);

		// calculate the fiber tangent
		c = NxN * (cn / J);
	}
	else
	{
		c.zero();
	}

	return c;
}

double FEFiberPow::FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
	double sed = 0.0;

	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// initialize material constants
	double ksi = m_ksi(mp);
	double beta = m_beta(mp);

	// loop over all integration points
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();

	// Calculate In = n0*C*n0
	double In = n0 * (C * n0);

	// only take fibers in tension into consideration
	if (!m_tension_only || (In - 1 >= eps))
	{
		// calculate strain energy density
		sed = ksi * pow(In - 1, beta);
	}

	return sed;
}

// define the material parameters
BEGIN_FECORE_CLASS(FEElasticFiberPow, FEElasticFiberMaterial)
	ADD_PARAMETER(m_fib.m_ksi, FE_RANGE_GREATER(0.0), "ksi")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_fib.m_beta, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETER(m_fib.m_tension_only, "tension_only");
END_FECORE_CLASS();
