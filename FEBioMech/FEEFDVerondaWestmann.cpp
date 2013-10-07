#include "stdafx.h"
#include "FEEFDVerondaWestmann.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEEFDVerondaWestmann, FEUncoupledMaterial)
	ADD_PARAMETER(m_VW.m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_VW.m_c2, FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETERV(m_EFD.m_beta, FE_PARAM_DOUBLEV, 3, "beta");
	ADD_PARAMETERV(m_EFD.m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEEFDVerondaWestmann::Init()
{
	FEUncoupledMaterial::Init();
	m_EFD.m_unstable = false;
	m_VW.Init();
	m_EFD.Init();
}

//-----------------------------------------------------------------------------
mat3ds FEEFDVerondaWestmann::DevStress(FEMaterialPoint& pt)
{
	return m_VW.DevStress(pt) + m_EFD.DevStress(pt);
}

//-----------------------------------------------------------------------------
tens4ds FEEFDVerondaWestmann::DevTangent(FEMaterialPoint& pt)
{
	return m_VW.DevTangent(pt) + m_EFD.DevTangent(pt);
}
