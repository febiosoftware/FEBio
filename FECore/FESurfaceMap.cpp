#include "stdafx.h"
#include "FESurfaceMap.h"
#include "FESurface.h"
#include "DumpStream.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap(int dataType) : FEDataArray(dataType)
{
	m_maxFaceNodes = 0;
}

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap(const FESurfaceMap& map) : FEDataArray(map), m_name(map.m_name)
{
	m_maxFaceNodes = map.m_maxFaceNodes;
}

//-----------------------------------------------------------------------------
FESurfaceMap& FESurfaceMap::operator = (const FESurfaceMap& map)
{
	FEDataArray::operator=(map);
	m_name = map.m_name;
	m_maxFaceNodes = map.m_maxFaceNodes;
	return *this;
}

//-----------------------------------------------------------------------------
bool FESurfaceMap::Create(const FESurface* ps, double val)
{
	int NF = ps->Elements();
	m_maxFaceNodes = 0;
	for (int i=0; i<NF; ++i)
	{
		const FESurfaceElement& el = ps->Element(i);
		int nf = el.Nodes();
		if (nf > m_maxFaceNodes) m_maxFaceNodes = nf;
	}
	return resize(NF*m_maxFaceNodes, val);
}

//-----------------------------------------------------------------------------
bool FESurfaceMap::Create(const FEFacetSet* ps, double val)
{
	int NF = ps->Faces();
	m_maxFaceNodes = 0;
	for (int i = 0; i<NF; ++i)
	{
		const FEFacetSet::FACET& f = ps->Face(i);

		// TODO: currently, the number of nodes matches the type, but not sure if this will remain the case.
		if (f.ntype > m_maxFaceNodes) m_maxFaceNodes = f.ntype;
	}
	return resize(NF*m_maxFaceNodes, val);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::SetName(const std::string& name)
{
	m_name = name;
}

//-----------------------------------------------------------------------------
void FESurfaceMap::setValue(int n, double v)
{
	int index = n*m_maxFaceNodes;
	for (int i=0; i<m_maxFaceNodes; ++i) set<double>(index+i, v);	
}

//-----------------------------------------------------------------------------
void FESurfaceMap::setValue(int n, const vec2d& v)
{
	int index = n*m_maxFaceNodes;
	for (int i = 0; i<m_maxFaceNodes; ++i) set<vec2d>(index + i, v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::setValue(int n, const vec3d& v)
{
	int index = n*m_maxFaceNodes;
	for (int i = 0; i<m_maxFaceNodes; ++i) set<vec3d>(index + i, v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::fillValue(double v)
{
	set<double>(v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::fillValue(const vec2d& v)
{
	set<vec2d>(v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::fillValue(const vec3d& v)
{
	set<vec3d>(v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::Serialize(DumpStream& ar)
{
	FEDataArray::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_maxFaceNodes;
		ar << m_name;
	}
	else
	{
		ar >> m_maxFaceNodes;
		ar >> m_name;
	}
}
