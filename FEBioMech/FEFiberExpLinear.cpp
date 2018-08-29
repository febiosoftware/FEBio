#include "stdafx.h"
#include "FEFiberExpLinear.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberExpLinear, FEElasticMaterial)
	ADD_PARAMETER2(m_c3  , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "c3");
	ADD_PARAMETER2(m_c4  , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "c4");
	ADD_PARAMETER2(m_c5  , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "c5");
	ADD_PARAMETER2(m_lam1, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(1.0), "lambda");
END_PARAMETER_LIST();

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
mat3ds FEFiberExpLinear::Stress(FEMaterialPoint& mp)
{
	// get the material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the material fiber axis
	vec3d a0 = GetFiberVector(mp);

	// get the spatial fiber vector and stretch
	vec3d a = pt.m_F*a0;
	double l = a.unit();

	// other stuff we need
	mat3ds A = dyad(a);
	double J = pt.m_J;

	// fiber stress
	mat3ds s; s.zero();

	// calculate fiber stress
	if (l > 1.0)
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
tens4ds FEFiberExpLinear::Tangent(FEMaterialPoint& mp)
{
	// get material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the material fiber axis
	vec3d a0 = GetFiberVector(mp);

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
	if (l > 1.0)
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
double FEFiberExpLinear::StrainEnergyDensity(FEMaterialPoint& mp)
{
	// TODO: implement this
	return 0;
}
