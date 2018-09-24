#include "stdafx.h"
#include "FETrussMaterial.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FETrussMaterial, FEMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// Note that this function returns the Kirchhoff stress!
double FETrussMaterial::Stress(FEMaterialPoint &mp)
{
	FETrussMaterialPoint& pt = *mp.ExtractData<FETrussMaterialPoint>();
	return m_E*log(pt.m_l);
}

//-----------------------------------------------------------------------------
double FETrussMaterial::Tangent(FEMaterialPoint &pt)
{
	return m_E;
}
