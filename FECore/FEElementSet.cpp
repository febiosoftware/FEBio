#include "stdafx.h"
#include "FEElementSet.h"
#include "FEElement.h"
#include "FEMesh.h"
#include "FEDomain.h"
#include "DumpStream.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
FEElementSet::FEElementSet(FEMesh* pm) : m_mesh(pm)
{
	m_minID = -1;
	m_maxID = -1;
}

//-----------------------------------------------------------------------------
void FEElementSet::Create(const std::vector<int>& elemList)
{
	m_dom.Clear();
	m_Elem = elemList;
	BuildLUT();
}

//-----------------------------------------------------------------------------
void FEElementSet::Create(FEDomain* dom)
{
	m_dom.Clear();
	m_dom.AddDomain(dom);

	int NE = dom->Elements();
	m_Elem.resize(NE, -1);
	for (int i = 0; i < NE; ++i)
	{
		FEElement& el = dom->ElementRef(i);
		m_Elem[i] = el.GetID();
	}

	BuildLUT();
}

//-----------------------------------------------------------------------------
void FEElementSet::Create(FEDomainList& domList)
{
	int NT = 0;
	m_dom.Clear();
	for (int n = 0; n < domList.Domains(); ++n)
	{
		FEDomain* dom = domList.GetDomain(n);
		m_dom.AddDomain(dom);
		NT += dom->Elements();
	}

	m_Elem.resize(NT, -1);
	NT = 0;
	for (int n = 0; n < domList.Domains(); ++n)
	{
		FEDomain* dom = domList.GetDomain(n);
		int NE = dom->Elements();
		for (int i = 0; i < NE; ++i)
		{
			FEElement& el = dom->ElementRef(i);
			m_Elem[NT + i] = el.GetID();
		}
		NT += NE;
	}

	BuildLUT();
}

//-----------------------------------------------------------------------------
void FEElementSet::BuildLUT()
{
	int N = (int)m_Elem.size();
	m_minID = m_maxID = -1;
	for (int i = 0; i < N; ++i)
	{
		FEElement* pe = m_mesh->FindElementFromID(m_Elem[i]);
		int id = pe->GetID();

		if ((id < m_minID) || (m_minID == -1)) m_minID = id;
		if ((id > m_maxID) || (m_maxID == -1)) m_maxID = id;
	}

	int lutSize = m_maxID - m_minID + 1;
	m_LUT.resize(lutSize, -1);
	for (int i = 0; i < N; ++i)
	{
		FEElement* pe = m_mesh->FindElementFromID(m_Elem[i]);
		int id = pe->GetID() - m_minID;
		m_LUT[id] = i;
	}
}

//-----------------------------------------------------------------------------
FEElement& FEElementSet::Element(int i)
{
	return *m_mesh->FindElementFromID(m_Elem[i]);
}

//-----------------------------------------------------------------------------
void FEElementSet::Serialize(DumpStream& ar)
{
	FEItemList::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_Elem;
	ar & m_LUT;
	ar & m_minID & m_maxID;
}

//-----------------------------------------------------------------------------
void FEElementSet::SaveClass(DumpStream& ar, FEElementSet* p)
{
}

//-----------------------------------------------------------------------------
FEElementSet* FEElementSet::LoadClass(DumpStream& ar, FEElementSet* p)
{
	FEMesh* mesh = &ar.GetFEModel().GetMesh();
	return new FEElementSet(mesh);
}
