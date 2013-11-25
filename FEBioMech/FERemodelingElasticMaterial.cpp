#include "stdafx.h"
#include "FERemodelingElasticMaterial.h"

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
	FEElasticMaterialPoint& pt = *m_pt->ExtractData<FEElasticMaterialPoint>();
	if (bflag)
	{
		// intialize data to zero
        dsed = rhorp = 0;
	}
	else
	{
		rhorp = pt.m_rhor;
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
		dmp << dsed << rhorp;
	}
	else
	{
		dmp >> dsed >> rhorp;
	}
}

//-----------------------------------------------------------------------------
void FERemodelingMaterialPoint::Serialize(DumpFile& ar)
{
	if (m_pt) m_pt->Serialize(ar);
        
	if (ar.IsSaving())
	{
		ar << dsed << rhorp;
	}
	else
	{
		ar >> dsed >> rhorp;
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
FEMaterial* FERemodelingElasticMaterial::GetProperty(int i)
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
bool FERemodelingElasticMaterial::SetProperty(int n, FEMaterial* pm)
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
//! Strain energy density function
double FERemodelingElasticMaterial::StrainEnergy(FEMaterialPoint& mp)
{
	return m_pBase->StrainEnergy(mp);
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FERemodelingElasticMaterial::Stress(FEMaterialPoint& mp)
{
	double dt = FEMaterialPoint::dt;

    FERemodelingMaterialPoint& rpt = *(mp.ExtractData<FERemodelingMaterialPoint>());
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

	// calculate the sed derivative with respect to mass density at this material point
    rpt.dsed = Tangent_SE_Density(mp);
                
	double rhorhat = m_pSupp->Supply(mp);
	pt.m_rhor = rhorhat*dt + rpt.rhorp;
	if (pt.m_rhor > m_rhormax) pt.m_rhor = m_rhormax;
	if (pt.m_rhor < m_rhormin) pt.m_rhor = m_rhormin;

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
