#include "stdafx.h"
#include "FEDomain.h"
#include "FEMaterial.h"
#include <string.h>

//-----------------------------------------------------------------------------
FEDomain::FEDomain(int ntype, int nclass, FEMesh* pm) : m_pMesh(pm), m_ntype(ntype), m_nclass(nclass)
{
	m_szname[0] = 0;
}

//-----------------------------------------------------------------------------
FEDomain::~FEDomain()
{
	
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
