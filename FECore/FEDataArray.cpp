#include "stdafx.h"
#include "FEDataArray.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
template <> bool FEDataArray::set<double>(int n, const double& v)
{
	assert(m_dataType == FE_DOUBLE);
	m_val[n] = v;
	return true;
}

//-----------------------------------------------------------------------------
template <> bool FEDataArray::set<vec2d>(int n, const vec2d& v)
{
	assert(m_dataType == FE_VEC2D);
	m_val[2*n  ] = v.x();
	m_val[2*n+1] = v.y();
	return true;
}

//-----------------------------------------------------------------------------
template <> bool FEDataArray::set<vec3d>(int n, const vec3d& v)
{
	assert(m_dataType == FE_VEC3D);
	m_val[3*n  ] = v.x;
	m_val[3*n+1] = v.y;
	m_val[3*n+2] = v.z;
	return true;
}

//-----------------------------------------------------------------------------
template <> bool FEDataArray::set<double>(const double& v)
{
	assert(m_dataType == FE_DOUBLE);
	m_defDouble = v;
	for (int i=0; i<(int) m_val.size(); ++i) m_val[i] = v;
	return true;
}

//-----------------------------------------------------------------------------
template <> bool FEDataArray::set<vec2d>(const vec2d& v)
{
	assert(m_dataType == FE_VEC2D);
	m_defVec2d = v;
	for (int i=0; i<(int) m_val.size(); i += 2)
	{
		m_val[i  ] = v.x();
		m_val[i+1] = v.y();
	}
	return true;
}

//-----------------------------------------------------------------------------
template <> bool FEDataArray::set<vec3d>(const vec3d& v)
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
FEDataArray::FEDataArray(int dataType) : m_dataType(dataType)
{
	m_defDouble = 0.0;
	m_defVec3d = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
FEDataArray::FEDataArray(const FEDataArray& map)
{
	m_dataType = map.m_dataType;
	m_val = map.m_val;
}

//-----------------------------------------------------------------------------
FEDataArray& FEDataArray::operator = (const FEDataArray& map)
{
	m_dataType = map.m_dataType;
	m_val = map.m_val;
	return *this;
}

//-----------------------------------------------------------------------------
int FEDataArray::DataSize() const
{
	switch (m_dataType)
	{
	case FE_DOUBLE : return 1; break;
	case FE_VEC2D  : return 2; break;
	case FE_VEC3D  : return 3; break;
	}
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
int FEDataArray::size() const
{
	return (int)m_val.size() / DataSize();
}

//-----------------------------------------------------------------------------
bool FEDataArray::resize(int n)
{
	switch (m_dataType)
	{
	case FE_DOUBLE: m_val.resize(n*DataSize(), m_defDouble); break;
	case FE_VEC2D : m_val.resize(n*DataSize()); set(m_defVec2d); break;
	case FE_VEC3D : m_val.resize(n*DataSize()); set(m_defVec3d); break;
	default:
		assert(false);
		return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
void FEDataArray::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_dataType;
		ar << m_val;
	}
	else
	{
		ar >> m_dataType;
		ar >> m_val;
	}
}
