#include "stdafx.h"
#include "FEElementSet.h"
#include "FEElement.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FEElementSet::FEElementSet(FEMesh* pm) : m_mesh(pm)
{
	m_minID = -1;
	m_maxID = -1;
}

//-----------------------------------------------------------------------------
void FEElementSet::Create(const std::vector<int>& elemList)
{
	m_Elem = elemList;
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
void FEElementSet::Serialize(DumpStream& ar)
{
	FEItemList::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_Elem;
		ar << m_LUT;
		ar << m_minID << m_maxID;
	}
	else
	{
		ar >> m_Elem;
		ar >> m_LUT;
		ar >> m_minID >> m_maxID;
	}
}
