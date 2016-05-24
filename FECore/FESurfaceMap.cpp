#include "stdafx.h"
#include "FESurfaceMap.h"
#include "FESurface.h"
#include "DumpStream.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
template <> bool FESurfaceMap::SetValue<double>(const FEFacetIndex& n, const double& v)
{
	assert(m_dataType == FE_DOUBLE);
	m_val[n] = v;
	return true;
}

//-----------------------------------------------------------------------------
template <> bool FESurfaceMap::SetValue<vec3d>(const FEFacetIndex& n, const vec3d& v)
{
	assert(m_dataType == FE_VEC3D);
	m_val[3*n  ] = v.x;
	m_val[3*n+1] = v.y;
	m_val[3*n+2] = v.z;
	return true;
}

//-----------------------------------------------------------------------------
template <> bool FESurfaceMap::SetValue<double>(const double& v)
{
	assert(m_dataType == FE_DOUBLE);
	m_defDouble = v;
	for (int i=0; i<(int) m_val.size(); ++i) m_val[i] = v;
	return true;
}

//-----------------------------------------------------------------------------
template <> bool FESurfaceMap::SetValue<vec3d>(const vec3d& v)
{
	assert(m_dataType == FE_VEC3D);
	m_defVec3d = v;
	for (int i=0; i<(int) m_val.size(); i += 3)
	{
		m_val[i  ] = v.x;
		m_val[i+1] = v.y;
		m_val[i+2] = v.z;
	}
	return true;
}

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap(int dataType) : m_dataType(dataType)
{
	m_defDouble = 0.0;
	m_defVec3d = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap(const FESurfaceMap& map)
{
	m_dataType = map.m_dataType;
	m_val = map.m_val;
	m_name = map.m_name;
}

//-----------------------------------------------------------------------------
FESurfaceMap& FESurfaceMap::operator = (const FESurfaceMap& map)
{
	m_dataType = map.m_dataType;
	m_val = map.m_val;
	m_name = map.m_name;
	return *this;
}

//-----------------------------------------------------------------------------
int FESurfaceMap::DataSize() const
{
	switch (m_dataType)
	{
	case FE_DOUBLE : return 1; break;
	case FE_VEC3D  : return 3; break;
	}
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
bool FESurfaceMap::Create(const FESurface* ps)
{
	int NF = ps->Elements();
	switch (m_dataType)
	{
	case FE_DOUBLE: m_val.resize(NF*DataSize(), m_defDouble); break;
	case FE_VEC3D : m_val.resize(NF*DataSize()); SetValue(m_defVec3d); break;
	default:
		assert(false);
		return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FESurfaceMap::Create(const FEFacetSet* ps)
{
	int NF = ps->Faces();
	switch (m_dataType)
	{
	case FE_DOUBLE: m_val.resize(NF*DataSize(), m_defDouble); break;
	case FE_VEC3D : m_val.resize(NF*DataSize()); SetValue(m_defVec3d); break;
	default:
		assert(false);
		return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
void FESurfaceMap::SetName(const std::string& name)
{
	m_name = name;
}

//-----------------------------------------------------------------------------
void FESurfaceMap::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_dataType;
		ar << m_name;
		ar << m_val;
	}
	else
	{
		ar >> m_dataType;
		ar >> m_name;
		ar >> m_val;
	}
}
