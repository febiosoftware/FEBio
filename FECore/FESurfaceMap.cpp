#include "stdafx.h"
#include "FESurfaceMap.h"
#include "FESurface.h"
#include "DumpStream.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap(int dataType) : FEDataArray(dataType)
{
}

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap(const FESurfaceMap& map) : FEDataArray(map), m_name(map.m_name)
{
}

//-----------------------------------------------------------------------------
FESurfaceMap& FESurfaceMap::operator = (const FESurfaceMap& map)
{
	FEDataArray::operator=(map);
	m_name = map.m_name;
	return *this;
}

//-----------------------------------------------------------------------------
bool FESurfaceMap::Create(const FESurface* ps)
{
	int NF = ps->Elements();
	return resize(NF);
}

//-----------------------------------------------------------------------------
bool FESurfaceMap::Create(const FEFacetSet* ps)
{
	int NF = ps->Faces();
	return resize(NF);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::SetName(const std::string& name)
{
	m_name = name;
}

//-----------------------------------------------------------------------------
void FESurfaceMap::Serialize(DumpStream& ar)
{
	FEDataArray::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_name;
	}
	else
	{
		ar >> m_name;
	}
}
