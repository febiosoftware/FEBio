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
#include "FEFiberNHUC.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberNHUC, FEFiberMaterialUncoupled)
	ADD_PARAMETER(m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu")->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFiberNHUC::FEFiberNHUC(FEModel* pfem) : FEFiberMaterialUncoupled(pfem) 
{ 
	m_mu = 0; 
}

//-----------------------------------------------------------------------------
mat3ds FEFiberNHUC::DevFiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J, -1.0 / 3.0);

	// loop over all integration points
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;

	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;

	// only take fibers in tension into consideration
	if (In_1 >= eps)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;

		// calculate the outer product of nt
		mat3ds N = dyad(nt);

		// calculate the fiber stress
		s = N*(m_mu*In_1 / J);
	}
	else
	{
		s.zero();
	}

	return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEFiberNHUC::DevFiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J, -1.0 / 3.0);

	// loop over all integration points
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;
	tens4ds c;

	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;

	// only take fibers in tension into consideration
	if (In_1 >= eps)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;

		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		tens4ds NxN = dyad1s(N);

		// calculate the fiber stress
		s = N*(m_mu*In_1 / J);

		// calculate the fiber tangent
		c = NxN*(2 * m_mu / J);

        // This is the final value of the elasticity tensor
        mat3dd I(1);
        tens4ds IxI = dyad1s(I);
        tens4ds I4 = dyad4s(I);
        c += ((I4 + IxI / 3.0)*s.tr() - dyad1s(I, s))*(2. / 3.)
        - (ddots(IxI, c) - IxI*(c.tr() / 3.)) / 3.;
	}
	else
	{
		c.zero();
	}

	return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberNHUC::DevFiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
	double sed = 0.0;

	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// loop over all integration points
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();

	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;

	// only take fibers in tension into consideration
	if (In_1 >= eps)
		sed = 0.25*m_mu*In_1*In_1;

	return sed;
}

// define the material parameters
BEGIN_FECORE_CLASS(FEUncoupledFiberNH, FEElasticFiberMaterialUC)
	ADD_PARAMETER(m_fib.m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu")->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();
