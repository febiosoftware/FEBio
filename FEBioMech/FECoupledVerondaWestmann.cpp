#include "stdafx.h"
#include "FECoupledVerondaWestmann.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FECoupledVerondaWestmann, FEElasticMaterial)
	ADD_PARAMETER(m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c2, FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETER(m_k , FE_PARAM_DOUBLE, "k" );
END_PARAMETER_LIST();


//-----------------------------------------------------------------------------
//! data initialization
void FECoupledVerondaWestmann::Init()
{
	FEElasticMaterial::Init();
	if (m_c1 <= 0.0) throw MaterialError("c1 must be positive");
	if (m_k  <= 0.0) throw MaterialError("k must be positive" );
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FECoupledVerondaWestmann::Stress(FEMaterialPoint& mp)
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
	return B*(m_c1*m_c2*(2.0*exp(m_c2*(I1-3)) - I1)/J) + B2*(m_c1*m_c2/J) + I*(m_k*log(J)/J);
}

//-----------------------------------------------------------------------------
//! calculate tangent at material point
tens4ds FECoupledVerondaWestmann::Tangent(FEMaterialPoint& mp)
{
	// get material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.LeftCauchyGreen();

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double J = pt.m_J;

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
	tens4ds c = BxB*(2.0*W1*W2*(2.0*W2*exp(W2*(I1-3))-1)/J) + BoB*(2.0*W1*W2/J) + IxI*(m_k/J) - IoI*(2.0*m_k*log(J)/J);

	return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FECoupledVerondaWestmann::StrainEnergyDensity(FEMaterialPoint& mp)
{
	// get the elastic material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// determinant of deformation gradient
	double J = pt.m_J;
    double lnJ = log(J);
    
	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();
    
	// calculate square of B
	mat3ds B2 = B*B;
    
	// Invariants of B (= invariants of C)
	double I1 = B.tr();
    double I2 = (I1*I1 - B2.tr())/2.0;
    
    double sed = m_c1*(exp(m_c2*(I1-3))-1) - m_c1*m_c2*(I2-3)/2 + m_k*lnJ*lnJ/2;
    
    return sed;
}

