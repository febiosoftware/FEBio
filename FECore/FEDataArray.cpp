#include "stdafx.h"
#include "FEDataArray.h"
#include "DumpStream.h"
#include "fecore_type.h"

//-----------------------------------------------------------------------------
FEDataArray::FEDataArray(FEDataMapType mapType, FEDataType dataType) : m_mapType(mapType), m_dataType(dataType), m_dataCount(0)
{
	m_dataSize = fecore_data_size(dataType);
}

//-----------------------------------------------------------------------------
FEDataArray::~FEDataArray()
{
}

//-----------------------------------------------------------------------------
FEDataArray::FEDataArray(const FEDataArray& map)
{
	m_dataType = map.m_dataType;
	m_dataSize = map.m_dataSize;
	m_dataCount = map.m_dataCount;
	m_val = map.m_val;
}

//-----------------------------------------------------------------------------
FEDataArray& FEDataArray::operator = (const FEDataArray& map)
{
	m_dataType = map.m_dataType;
	m_dataSize = map.m_dataSize;
	m_dataCount = map.m_dataCount;
	m_val = map.m_val;
	return *this;
}

//-----------------------------------------------------------------------------
bool FEDataArray::resize(int n, double val)
{
	if (n < 0) return false;
	m_val.resize(n*DataSize(), val);
	m_dataCount = n;
	return true;
}

//-----------------------------------------------------------------------------
//! set the data sized
void FEDataArray::SetDataSize(int dataSize)
{
	m_dataSize = dataSize;
	if (m_val.empty() == false)
	{
		m_val.resize(m_dataSize*m_dataCount);
	}
}

//-----------------------------------------------------------------------------

FEDataType intToDataType(int i)
{
	switch (i)
	{
	case FE_DOUBLE: return FE_DOUBLE; break;
	case FE_VEC2D : return FE_VEC2D; break;
	case FE_VEC3D : return FE_VEC3D; break;
	case FE_MAT3D : return FE_MAT3D; break;
	default:
		assert(false);
		break;
	}
	return FE_INVALID_TYPE;
}

void FEDataArray::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_dataSize;
		ar << m_dataCount;
		ar << (int)m_dataType;
		ar << m_val;
	}
	else
	{
		int ntype = 0;
		ar >> m_dataSize;
		ar >> m_dataCount;
		ar >> ntype;
		ar >> m_val;

		m_dataType = intToDataType(ntype);
		assert(m_val.size() == m_dataSize*m_dataCount);
	}
}
