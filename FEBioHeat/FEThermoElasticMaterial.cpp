#include "FEThermoElasticMaterial.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//! FEBiphasic constructor

FEThermoElasticMaterial::FEThermoElasticMaterial(FEModel* pfem) : FEMaterial(pfem)
{ 
	// set material properties
	m_pElastic.SetName("elastic").SetID(0);
	m_pCond.SetName("conductivity").SetID(1);
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEThermoElasticMaterial::CreateMaterialPointData() 
{ 
	return new FEHeatMaterialPoint(m_pElastic->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
void FEThermoElasticMaterial::Init()
{
	FEMaterial::Init();
	m_pElastic->SetParent(this); m_pElastic->Init();
	m_pCond->SetParent(this); m_pCond->Init();
}

//-----------------------------------------------------------------------------
int FEThermoElasticMaterial::MaterialProperties()
{
	return 2;
}

//-----------------------------------------------------------------------------
//! return a pointer to a biphasic material property
FEProperty* FEThermoElasticMaterial::GetMaterialProperty(int i)
{
	switch (i)
	{
	case 0: return &m_pElastic;
	case 1: return &m_pCond;
	}
	assert(false);
	return 0;
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
