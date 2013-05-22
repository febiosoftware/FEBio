#include "stdafx.h"
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
// Material parameters for FEUncoupledMaterial
BEGIN_PARAMETER_LIST(FEUncoupledMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_K, FE_PARAM_DOUBLE, "k");
	ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL  , "laugon");
	ADD_PARAMETER(m_atol   , FE_PARAM_DOUBLE, "atol"  );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEUncoupledMaterial::FEUncoupledMaterial(void)
{
	m_blaugon = false;
	m_atol = 0.01;
	m_K = 0;	// invalid value!
}

//-----------------------------------------------------------------------------
void FEUncoupledMaterial::Init()
{
	FEElasticMaterial::Init();
	if (m_K < 0) throw MaterialError("k must be positive.");
}

//-----------------------------------------------------------------------------
//! The stress function calculates the total Cauchy stress as a sum of 
//! two terms, namely the deviatoric stress and the pressure. 
mat3ds FEUncoupledMaterial::Stress(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate the stress as a sum of deviatoric stress and pressure
	return mat3dd(UJ(pt.J)) + DevStress(mp);
}

//------------------------------------------------------------------------------
//! The tangent function calculates the total spatial tangent, that is it calculates
//! the push-forward of the derivative of the 2ndPK stress with respect to C. However,
//! for an uncoupled material, the 2ndPK stress decouples in a deviatoric and a 
//! dilatational component. The deviatoric tangent is provided by the particular
//! material and the dilatational component is added here.
//!
tens4ds FEUncoupledMaterial::Tangent(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// 2nd-order identity tensor
	mat3dd I(1);

	// 4th-order identity tensors
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	
	// pressure
	double p = UJ(pt.J);
	
	// tangent is sum of three terms
	// C = c_tilde + c_pressure + c_k
	//
	// + c_tilde is the derivative of the deviatoric stress with respect to C
	// + c_pressure is p*d(JC)/dC
	// + c_k comes from the derivative of p with respect to C
	// 
	// Note that the c_k term is not necessary in the 3F formulation (since p is independant variable) 
	// but we do need to add it here.
	//
	//        c_tilde         c_pressure            c_k
	return DevTangent(mp) + (IxI - I4*2)*p + IxI*(UJJ(pt.J)*pt.J);
}
