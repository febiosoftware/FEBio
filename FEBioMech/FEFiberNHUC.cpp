#include "stdafx.h"
#include "FEFiberNHUC.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberNHUC, FEElasticFiberMaterialUC)
	ADD_PARAMETER(m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
mat3ds FEFiberNHUC::DevStress(FEMaterialPoint& mp, const vec3d& n0)
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
	if (In_1 > eps)
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
tens4ds FEFiberNHUC::DevTangent(FEMaterialPoint& mp, const vec3d& n0)
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
	if (In_1 > eps)
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
	}
	else
	{
		c.zero();
	}

	// This is the final value of the elasticity tensor
	mat3dd I(1);
	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);
	c += ((I4 + IxI / 3.0)*s.tr() - dyad1s(I, s))*(2. / 3.)
		- (ddots(IxI, c) - IxI*(c.tr() / 3.)) / 3.;

	return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberNHUC::DevStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
	double sed = 0.0;

	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// loop over all integration points
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();

	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;

	// only take fibers in tension into consideration
	if (In_1 > eps)
		sed = 0.25*m_mu*In_1*In_1;

	return sed;
}
