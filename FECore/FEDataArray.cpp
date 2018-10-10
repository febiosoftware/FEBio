#include "stdafx.h"
#include "FEDataArray.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
FEDataArray::FEDataArray(int dataType) : m_dataSize(dataType), m_dataCount(0)
{
}

//-----------------------------------------------------------------------------
FEDataArray::~FEDataArray()
{
}

//-----------------------------------------------------------------------------
FEDataArray::FEDataArray(const FEDataArray& map)
{
	m_dataSize = map.m_dataSize;
	m_dataCount = map.m_dataCount;
	m_val = map.m_val;
}

//-----------------------------------------------------------------------------
FEDataArray& FEDataArray::operator = (const FEDataArray& map)
{
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
void FEDataArray::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_dataSize;
		ar << m_dataCount;
		ar << m_val;
	}
	else
	{
		ar >> m_dataSize;
		ar >> m_dataCount;
		ar >> m_val;
		assert(m_val.size() == m_dataSize*m_dataCount);
	}
}
