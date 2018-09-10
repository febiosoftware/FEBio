#include "stdafx.h"
#include "FESegmentSet.h"

//-----------------------------------------------------------------------------
FESegmentSet::FESegmentSet(FEMesh* pm) : m_mesh(pm)
{
}

//-----------------------------------------------------------------------------
void FESegmentSet::Create(int n)
{
	m_Seg.resize(n);
}

//-----------------------------------------------------------------------------
FESegmentSet::SEGMENT& FESegmentSet::Segment(int i)
{
	return m_Seg[i];
}

//-----------------------------------------------------------------------------
void FESegmentSet::SetName(const std::string& name)
{
	m_name = name;
}

//-----------------------------------------------------------------------------
const std::string& FESegmentSet::GetName() const
{
	return m_name;
}

//-----------------------------------------------------------------------------
void FESegmentSet::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_name;
		ar << m_Seg;
	}
	else
	{
		ar >> m_name;
		ar >> m_Seg;
	}
}
