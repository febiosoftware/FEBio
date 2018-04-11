#include "stdafx.h"
#include "FEDiscreteDomain.h"
#include "FEMesh.h"
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
void FEDiscreteDomain::Create(int nelems, int elemType)
{ 
	m_Elem.resize(nelems); 
	for (int i = 0; i<nelems; ++i) m_Elem[i].SetDomain(this);

	if (elemType != -1)
		for (int i=0; i<nelems; ++i) m_Elem[i].SetType(elemType);
}

//-----------------------------------------------------------------------------
bool FEDiscreteDomain::Init()
{
	if (FEDomain::Init() == false) return false;

	FEMaterial* pmat = GetMaterial();
	if (pmat) SetMatID(pmat->GetID());

	return true;
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::AddElement(int eid, int n[2])
{
	FEDiscreteElement el;
	el.SetType(FE_DISCRETE);
	el.m_node[0] = n[0];
	el.m_node[1] = n[1];
	el.SetID(eid);
	m_Elem.push_back(el);
}
