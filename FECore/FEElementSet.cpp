#include "stdafx.h"
#include "FEElementSet.h"

//-----------------------------------------------------------------------------
FEElementSet::FEElementSet(FEMesh* pm) : m_mesh(pm)
{
}

//-----------------------------------------------------------------------------
void FEElementSet::create(int n)
{
	assert(n);
	m_Elem.resize(n);
}

//-----------------------------------------------------------------------------
void FEElementSet::SetName(const std::string& name)
{
	m_name = name;
}

//-----------------------------------------------------------------------------
const std::string& FEElementSet::GetName() const
{
	return m_name;
}

//-----------------------------------------------------------------------------
void FEElementSet::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_name;
		ar << m_Elem;
	}
	else
	{
		ar >> m_name;
		ar >> m_Elem;
	}
}
