/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEDataArray.h"
#include "DumpStream.h"
#include "fecore_type.h"
#include "FENodeDataMap.h"
#include "FEDomainMap.h"
#include "FESurfaceMap.h"

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
bool FEDataArray::realloc(int n)
{
	if (n < 0) return false;
	m_val.resize(n*DataSize());
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
    case FE_MAT3DS: return FE_MAT3DS; break;
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

void FEDataArray::SaveClass(DumpStream& ar, FEDataArray* p)
{
	int ntype = (int) p->DataMapType();
	ar << ntype;
}

FEDataArray* FEDataArray::LoadClass(DumpStream& ar, FEDataArray* p)
{
	int ntype;
	ar >> ntype;
	p = nullptr;
	switch (ntype)
	{
	case FE_NODE_DATA_MAP: p = new FENodeDataMap; break;
	case FE_DOMAIN_MAP   : p = new FEDomainMap; break;
	case FE_SURFACE_MAP  : p = new FESurfaceMap; break;
	default:
		assert(false);
	}

	return p;
}
