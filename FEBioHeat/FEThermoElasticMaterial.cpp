#include "FEThermoElasticMaterial.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//! FEBiphasic constructor

FEThermoElasticMaterial::FEThermoElasticMaterial(FEModel* pfem) : FEMaterial(pfem)
{ 
	m_pElastic = 0;
	m_pCond    = 0;
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
int FEThermoElasticMaterial::Properties()
{
    int np = 2;
	return np;
}

//-----------------------------------------------------------------------------
//! return a pointer to a biphasic material property
FECoreBase* FEThermoElasticMaterial::GetProperty(int i)
{
	switch (i)
	{
	case 0: return m_pElastic;
	case 1: return m_pCond;
	}
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
//! Find the index of a material property
int FEThermoElasticMaterial::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "elastic"     ) == 0) return 0;
	if (strcmp(szname, "conductivity") == 0) return 1;
	return -1;
}

//-----------------------------------------------------------------------------
//! Set a material property
bool FEThermoElasticMaterial::SetProperty(int n, FECoreBase* pm)
{
	switch(n)
	{
	case 0:
		{
			FEThermalElastic* pme = dynamic_cast<FEThermalElastic*>(pm);
			if (pme) { m_pElastic = pme; return true; }
		}
		break;
	case 1: 
		{
			FEThermalConductivity* pmc = dynamic_cast<FEThermalConductivity*>(pm);
			if (pmc) { m_pCond = pmc; return true; }
		}
		break;
	}
	return false;
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

//-----------------------------------------------------------------------------
//! serialization
void FEThermoElasticMaterial::Serialize(DumpFile &ar)
{
	// serialize material parameters
	FEMaterial::Serialize(ar);

	// serialize sub-materials
	if (ar.IsSaving())
	{
		ar << m_pElastic->GetTypeStr();
		m_pElastic->Serialize(ar);

		ar << m_pCond->GetTypeStr();
		m_pCond->Serialize(ar);
	}
	else
	{
		char sz[256] = {0};

		ar >> sz;
		m_pElastic = dynamic_cast<FEThermalElastic*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
		assert(m_pElastic);
		m_pElastic->Serialize(ar);
		m_pElastic->Init();

		ar >> sz;
		m_pCond = dynamic_cast<FEThermalConductivity*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
		assert(m_pCond);
		m_pCond->Serialize(ar);
		m_pCond->Init();
	}
}

//-----------------------------------------------------------------------------
FEParam* FEThermoElasticMaterial::GetParameter(const ParamString& s)
{
	// see if this is a composite parameter
	if (s.count() == 1) return FEMaterial::GetParameter(s);

	// else find the component's parameter
	if      (s == "elastic"     ) return m_pElastic->GetParameter(s.next());
	else if (s == "conductivity") return m_pCond   ->GetParameter(s.next());
	else return 0;
}
