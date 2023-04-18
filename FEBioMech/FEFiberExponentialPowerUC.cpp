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
#include "FEFiberExponentialPowerUC.h"
#include <limits>
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberExponentialPowerUC, FEElasticFiberMaterialUC)
	ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
	ADD_PARAMETER(m_beta, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETER(m_ksi, "ksi");
	ADD_PARAMETER(m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFiberExponentialPowerUC::FEFiberExponentialPowerUC(FEModel* pfem) : FEElasticFiberMaterialUC(pfem)
{
	m_alpha = 0;
	m_beta = 2;
	m_ksi = 0;
	m_mu = 0;
    m_epsf = 0.0;
}

//-----------------------------------------------------------------------------
bool FEFiberExponentialPowerUC::Validate()
{
	if ((4 * m_ksi + 2 * m_mu) < 0) { feLogError("4*ksi+2*mu must be positive."); return false;	}
	return FEElasticFiberMaterialUC::Validate();
}

//-----------------------------------------------------------------------------
mat3ds FEFiberExponentialPowerUC::DevFiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J, -1.0 / 3.0);

	// loop over all integration points
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;

	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;

	// only take fibers in tension into consideration
    const double eps = m_epsf* std::numeric_limits<double>::epsilon();
	if (In_1 >= eps)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;

		// calculate the outer product of nt
		mat3ds N = dyad(nt);

		// calculate strain energy derivative
		double Wl = m_ksi*pow(In_1, m_beta - 1.0)*exp(m_alpha*pow(In_1, m_beta));

		// calculate the fiber stress
		s = N*(2.0*Wl / J);

		// add the contribution from shear
		mat3ds BmI = pt.DevLeftCauchyGreen() - mat3dd(1);
		s += (N*BmI).sym()*(m_mu / J);
	}
	else
	{
		s.zero();
	}

	return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEFiberExponentialPowerUC::DevFiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J, -1.0 / 3.0);

    const double eps = m_epsf*std::numeric_limits<double>::epsilon();
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
        mat3dd I(1);
        tens4ds IxI = dyad1s(I);
        tens4ds I4 = dyad4s(I);

		// calculate strain energy derivative
		double Wl = m_ksi*pow(In_1, m_beta - 1.0)*exp(m_alpha*pow(In_1, m_beta));

		// calculate the fiber stress
		s = N*(2.0*Wl / J);

        // add the contribution from shear
        mat3ds B = pt.DevLeftCauchyGreen();
        mat3ds BmI = B - I;
        s += (N*BmI).sym()*(m_mu / J);

		// calculate strain energy 2nd derivative
		double tmp = m_alpha*pow(In_1, m_beta);
		double Wll = m_ksi*pow(In_1, m_beta - 2.0)*((tmp + 1)*m_beta - 1.0)*exp(tmp);

		// calculate the fiber tangent
		c = NxN*(4.0*Wll / J);

		// add the contribution from shear
		c += dyad4s(N, B)*(m_mu / J);

        // This is the final value of the elasticity tensor
        c += ((I4 + IxI / 3.0)*s.tr() - dyad1s(I, s))*(2. / 3.) - (ddots(IxI, c) - IxI*(c.tr() / 3.)) / 3.;
	}
	else
	{
		c.zero();
	}

	return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberExponentialPowerUC::DevFiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// loop over all integration points
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds C2 = C.sqr();

	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;

	// only take fibers in tension into consideration
    const double eps = m_epsf*std::numeric_limits<double>::epsilon();
	double sed = 0.0;
	if (In_1 >= eps)
	{
		// calculate strain energy density
		if (m_alpha > 0)
			sed = m_ksi / (m_alpha*m_beta)*(exp(m_alpha*pow(In_1, m_beta)) - 1);
		else
			sed = m_ksi / m_beta*pow(In_1, m_beta);

		// add the contribution from shear
		sed += m_mu*(n0*(C2*n0) - 2 * In_1 - 1) / 4.0;
	}

	return sed;
}
