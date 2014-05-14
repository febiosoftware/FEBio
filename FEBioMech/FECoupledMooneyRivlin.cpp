#include "stdafx.h"
#include "FECoupledMooneyRivlin.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FECoupledMooneyRivlin, FEElasticMaterial)
	ADD_PARAMETER(m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c2, FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETER(m_K , FE_PARAM_DOUBLE, "k" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! data initialization
void FECoupledMooneyRivlin::Init()
{
	if (m_c1 <= 0.0) throw MaterialError("c1 must be positive");
	if (m_K  <= 0.0) throw MaterialError("k must be positive" );
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FECoupledMooneyRivlin::Stress(FEMaterialPoint& mp)
{
	// get the elastic material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;

	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// Invariants of B (= invariants of C)
	double I1 = B.tr();

	// identity tensor
	mat3ds I(1.0);

	// calculate stress
	return (B*(m_c1+I1*m_c2) - B2*m_c2 - I*(m_c1+2.0*m_c2))*(2.0/J) + I*(m_K*log(J)/J);
}

//-----------------------------------------------------------------------------
//! calculate tangent at material point
tens4ds FECoupledMooneyRivlin::Tangent(FEMaterialPoint& mp)
{
	// get material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.LeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double J = pt.m_J;
	double Ji = 1.0/J;

	// some useful tensors
	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds IoI = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds BoB = dyad4s(B);

	// strain energy derivates
	double W1 = m_c1;
	double W2 = m_c2;

	// spatial tangent
	tens4ds c = BxB*(4.0*W2/J) - BoB*(4.0*W2/J) + IoI*(4.0*(m_c1+2.0*m_c2)/J) + IxI*(m_K/J) - IoI*(2.0*m_K*log(J)/J);

	return c;
}
