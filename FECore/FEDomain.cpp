#include "stdafx.h"
#include "FEDomain.h"
#include "FEMaterial.h"
#include <string.h>

//-----------------------------------------------------------------------------
void FEDataExport::Serialize(vector<float>& d)
{
	if ((m_type == PLT_VEC3F)&&(m_fmt == FMT_NODE))
	{
		vector<vec3d>& v = *(static_cast<vector<vec3d>*>(m_pd));

		int n = (int) v.size();
		for (int i=0; i<n; ++i)
		{
			d.push_back((float) v[i].x);
			d.push_back((float) v[i].y);
			d.push_back((float) v[i].z);
		}
	}
}

//-----------------------------------------------------------------------------
FEDomain::FEDomain(int nclass, FEMesh* pm) : FECoreBase(FEDOMAIN_ID), m_pMesh(pm), m_nclass(nclass)
{
	m_szname[0] = 0;
}

//-----------------------------------------------------------------------------
FEDomain::~FEDomain()
{
	// delete all data export classes
	if (m_Data.empty() == false)
	{
		size_t ND = m_Data.size();
		for (size_t i=0; i<ND; ++i) delete m_Data[i];
		m_Data.clear();
	}
}

//-----------------------------------------------------------------------------
void FEDomain::AddDataExport(FEDataExport* pd)
{
	if (pd) m_Data.push_back(pd);
}

//-----------------------------------------------------------------------------
void FEDomain::SetName(const char* szname)
{
	if (szname == 0) return;
	int l = strlen(szname);
	if (l>=MAX_DOMAIN_NAME) l = MAX_DOMAIN_NAME - 1;
	if (l>0) strncpy(m_szname, szname, l);
	m_szname[l] = 0;
}

//-----------------------------------------------------------------------------
const char* FEDomain::GetName()
{
	return m_szname;
}

//-----------------------------------------------------------------------------
FEElement* FEDomain::FindElementFromID(int nid)
{
	for (int i=0; i<Elements(); ++i)
	{
		FEElement& el = ElementRef(i);
		if (el.m_nID == nid) return &el;
	}

	return 0;
}

//-----------------------------------------------------------------------------
void FEDomain::InitMaterialPointData()
{
	FEMaterial* pmat = GetMaterial();

	for (int i=0; i<Elements(); ++i)
	{
		FEElement& el = ElementRef(i);
		for (int k=0; k<el.GaussPoints(); ++k) el.SetMaterialPointData(pmat->CreateMaterialPointData(), k);
	}
}

//-----------------------------------------------------------------------------
void FEDomain::SetMatID(int mid)
{
	for (int i=0; i<Elements(); ++i) ElementRef(i).SetMatID(mid);
}
