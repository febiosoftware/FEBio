#include "stdafx.h"
#include "FERemodelingElasticMaterial.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FERemodelingElasticMaterial, FEElasticMaterial)
ADD_PARAMETER(m_rhormin, FE_PARAM_DOUBLE, "min_density");
ADD_PARAMETER(m_rhormax, FE_PARAM_DOUBLE, "max_density");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Initialization
void FERemodelingElasticMaterial::Init()
{
	FEElasticMaterial::Init();
	m_pBase->Init();
	m_pSupp->Init();
}

//-----------------------------------------------------------------------------
//! Strain energy density function
double FERemodelingElasticMaterial::StrainEnergy(FEMaterialPoint& mp)
{
	return m_pBase->StrainEnergy(mp);
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FERemodelingElasticMaterial::Stress(FEMaterialPoint& mp)
{
	return m_pBase->Stress(mp);
}

//-----------------------------------------------------------------------------
//! Tangent of stress with strain
tens4ds FERemodelingElasticMaterial::Tangent(FEMaterialPoint& mp)
{
	return m_pBase->Tangent(mp);
}

//-----------------------------------------------------------------------------
//! Tangent of strain energy density with mass density
double FERemodelingElasticMaterial::Tangent_SE_Density(FEMaterialPoint& pt)
{
    return m_pBase->Tangent_SE_Density(pt);
}

//-----------------------------------------------------------------------------
//! Tangent of stress with mass density
mat3ds FERemodelingElasticMaterial::Tangent_Stress_Density(FEMaterialPoint& pt)
{
    return m_pBase->Tangent_Stress_Density(pt);
}
