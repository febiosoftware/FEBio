// FEIncompNeoHookean.cpp: implementation of the FEIncompNeoHookean class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEIncompNeoHookean.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEIncompNeoHookean, FEUncoupledMaterial)
	ADD_PARAMETER(m_G, FE_PARAM_DOUBLE, "G");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEIncompNeoHookean::Init()
{
	FEUncoupledMaterial::Init();
	if (m_G <= 0) throw MaterialError("G must be positive.");
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric stress
mat3ds FEIncompNeoHookean::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;

	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();

	// calculate deviatoric stress
	return (B - mat3dd(B.tr()/3.))*(m_G/pow(J, 5.0/3.0));
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4ds FEIncompNeoHookean::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;

	// left cauchy-green matrix (i.e. the 'b' matrix)
	mat3ds B = pt.LeftCauchyGreen();

	// trace of b
	double Ib = B.tr();

	double muJ = m_G/pow(J, 5.0/3.0);

	mat3ds I(1,1,1,0,0,0);	// Identity

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxI = dyad1s(B, I); // = BxI + IxB

	return (I4*Ib -BxI + IxI*(Ib/3))*(2.0*muJ/3.0);
}
