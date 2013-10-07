#include "stdafx.h"
#include "FEEFDMooneyRivlin.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEEFDMooneyRivlin, FEUncoupledMaterial)
	ADD_PARAMETER(m_MR.c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_MR.c2, FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETERV(m_EFD.m_beta, FE_PARAM_DOUBLEV, 3, "beta");
	ADD_PARAMETERV(m_EFD.m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEEFDMooneyRivlin::FEEFDMooneyRivlin()
{

}

//-----------------------------------------------------------------------------
void FEEFDMooneyRivlin::Init()
{
	FEUncoupledMaterial::Init();
	m_EFD.m_unstable = false;
	m_MR.Init();
	m_EFD.Init();
}

//-----------------------------------------------------------------------------
mat3ds FEEFDMooneyRivlin::DevStress(FEMaterialPoint& pt)
{
	return m_MR.DevStress(pt) + m_EFD.DevStress(pt);
}

//-----------------------------------------------------------------------------
tens4ds FEEFDMooneyRivlin::DevTangent(FEMaterialPoint& pt)
{
	return m_MR.DevTangent(pt) + m_EFD.DevTangent(pt);
}
