#include "stdafx.h"
#include "FERemodelingElasticMaterial.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
FEMaterialPoint* FERemodelingMaterialPoint::Copy()
{
	FERemodelingMaterialPoint* pt = new FERemodelingMaterialPoint(*this);
	if (m_pt) pt->m_pt = m_pt->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FERemodelingMaterialPoint::Init(bool bflag)
{
	if (bflag)
	{
		// intialize data to zero
        m_sed = m_dsed = 0; 
		m_rhor = m_rhorp = 0;
	}
	else
	{
		m_rhorp = m_rhor;
	}
        
	// don't forget to intialize the nested data
	if (m_pt) m_pt->Init(bflag);
}

//-----------------------------------------------------------------------------
void FERemodelingMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (m_pt) m_pt->ShallowCopy(dmp, bsave);
        
	if (bsave)
	{
		dmp << m_sed << m_dsed;
		dmp << m_rhor << m_rhorp;
	}
	else
	{
		dmp >> m_sed >> m_dsed;
		dmp >> m_rhor >> m_rhorp;
	}
}

//-----------------------------------------------------------------------------
void FERemodelingMaterialPoint::Serialize(DumpFile& ar)
{
	if (m_pt) m_pt->Serialize(ar);
        
	if (ar.IsSaving())
	{
		ar << m_sed << m_dsed;
		ar << m_rhor << m_rhorp;
	}
	else
	{
		ar >> m_sed >> m_dsed;
		ar >> m_rhor >> m_rhorp;
	}
}

//=============================================================================
// FERemodelingElasticMaterial
//=============================================================================

//-----------------------------------------------------------------------------
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
//! This material has two properties
int FERemodelingElasticMaterial::Properties()
{
	return 2;
}

//-----------------------------------------------------------------------------
FECoreBase* FERemodelingElasticMaterial::GetProperty(int i)
{
	switch (i)
	{
	case 0: return m_pBase;
	case 1: return m_pSupp;
	}
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
//! Find the index of a material property
int FERemodelingElasticMaterial::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "solid" ) == 0) return 0;
	if (strcmp(szname, "supply") == 0) return 1;
	return -1;
}

//-----------------------------------------------------------------------------
//! Set a material property
bool FERemodelingElasticMaterial::SetProperty(int n, FECoreBase* pm)
{
	switch(n)
	{
	case 0:
		{
			FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pm);
			if (pme) { m_pBase = pme; return true; }
		}
		break;
	case 1: 
		{
			FESolidSupply* pms = dynamic_cast<FESolidSupply*>(pm);
			if (pms) { m_pSupp = pms; return true; }
		}
		break;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Serialization
void FERemodelingElasticMaterial::Serialize(DumpFile &ar)
{
	// serialize material parameters
	FEElasticMaterial::Serialize(ar);

	// serialize sub-materials
	if (ar.IsSaving())
	{
		ar << m_pBase->GetTypeStr();
		m_pBase->Serialize(ar);

		ar << m_pSupp->GetTypeStr();
		m_pSupp->Serialize(ar);
	}
	else
	{
		char sz[256] = {0};

		ar >> sz;
		m_pBase = dynamic_cast<FEElasticMaterial*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
		assert(m_pBase);
		m_pBase->Serialize(ar);
		m_pBase->Init();

		ar >> sz;
		m_pSupp = dynamic_cast<FESolidSupply*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
		assert (m_pSupp);
		m_pSupp->Serialize(ar);
		m_pSupp->Init();
	}
}

//-----------------------------------------------------------------------------
//! Strain energy density function
double FERemodelingElasticMaterial::StrainEnergyDensity(FEMaterialPoint& mp)
{
	return (dynamic_cast<FERemodelingInterface*>(m_pBase))->StrainEnergy(mp);
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FERemodelingElasticMaterial::Stress(FEMaterialPoint& mp)
{
	double dt = FEMaterialPoint::dt;

    FERemodelingMaterialPoint& rpt = *(mp.ExtractData<FERemodelingMaterialPoint>());

	// calculate the strain energy density at this material point
	rpt.m_sed = StrainEnergyDensity(mp);

	// calculate the sed derivative with respect to mass density at this material point
    rpt.m_dsed = Tangent_SE_Density(mp);
                
	double rhorhat = m_pSupp->Supply(mp);
	rpt.m_rhor = rhorhat*dt + rpt.m_rhorp;
	if (rpt.m_rhor > m_rhormax) rpt.m_rhor = m_rhormax;
	if (rpt.m_rhor < m_rhormin) rpt.m_rhor = m_rhormin;

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
    return (dynamic_cast<FERemodelingInterface*>(m_pBase))->Tangent_SE_Density(pt);
}

//-----------------------------------------------------------------------------
//! Tangent of stress with mass density
mat3ds FERemodelingElasticMaterial::Tangent_Stress_Density(FEMaterialPoint& pt)
{
    return (dynamic_cast<FERemodelingInterface*>(m_pBase))->Tangent_Stress_Density(pt);
}
