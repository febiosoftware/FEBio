#include "FEThermoElasticMaterial.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//! FEBiphasic constructor

FEThermoElasticMaterial::FEThermoElasticMaterial(FEModel* pfem) : FEMaterial(pfem)
{ 
	// set material properties
	AddProperty(&m_pElastic, "elastic"     );
	AddProperty(&m_pCond   , "conductivity");
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEThermoElasticMaterial::CreateMaterialPointData() 
{ 
	return new FEHeatMaterialPoint(m_pElastic->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
mat3ds FEThermoElasticMaterial::Conductivity(FEMaterialPoint& mp)
{
	return m_pCond->Conductivity(mp);
}

//-----------------------------------------------------------------------------
mat3ds FEThermoElasticMaterial::Stress(FEMaterialPoint& mp)
{
	return m_pElastic->Stress(mp);
}

//-----------------------------------------------------------------------------
tens4ds FEThermoElasticMaterial::Tangent(FEMaterialPoint& mp)
{
	return m_pElastic->Tangent(mp);
}

//-----------------------------------------------------------------------------
mat3ds FEThermoElasticMaterial::ThermalTangent(FEMaterialPoint& mp)
{
	return m_pElastic->ThermalTangent(mp);
}

//-----------------------------------------------------------------------------
tens4ds FEThermoElasticMaterial::ConductivityGradient(FEMaterialPoint& mp)
{
	return m_pCond->Tangent_Conductivity_Strain(mp);
}
