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
#include "FENodeDataMap.h"
#include "FENodeSet.h"
#include "FEMaterialPoint.h"

FENodeDataMap::FENodeDataMap() : FEDataMap(FE_NODE_DATA_MAP, FE_INVALID_TYPE)
{
	m_nodeSet = nullptr;
}

FENodeDataMap::FENodeDataMap(FEDataType dataType) : FEDataMap(FE_NODE_DATA_MAP, dataType)
{
	m_nodeSet = nullptr;
}

void FENodeDataMap::Create(const FENodeSet* nodeSet, double val)
{
	m_nodeSet = nodeSet;
	int nsize = nodeSet->Size();
	resize(nsize, val);
}

const FENodeSet* FENodeDataMap::GetNodeSet() const
{ 
	return m_nodeSet; 
}

double FENodeDataMap::getValue(int n) const
{
	return get<double>(n);
}

void FENodeDataMap::setValue(int n, double v)
{
	set<double>(n, v);
}

void FENodeDataMap::setValue(int n, const vec2d& v)
{
	set<vec2d>(n, v);
}

void FENodeDataMap::setValue(int n, const vec3d& v)
{
	set<vec3d>(n, v);
}

void FENodeDataMap::setValue(int n, const mat3d& v)
{
	set<mat3d>(n, v);
}

void FENodeDataMap::setValue(int n, const mat3ds& v)
{
	set<mat3ds>(n, v);
}

void FENodeDataMap::fillValue(double v)
{
	set<double>(v);
}

void FENodeDataMap::fillValue(const vec2d& v)
{
	set<vec2d>(v);
}

void FENodeDataMap::fillValue(const vec3d& v)
{
	set<vec3d>(v);
}

void FENodeDataMap::fillValue(const mat3d& v)
{
	set<mat3d>(v);
}

void FENodeDataMap::fillValue(const mat3ds& v)
{
	set<mat3ds>(v);
}

double FENodeDataMap::value(const FEMaterialPoint& mp)
{
	assert(mp.m_elem == nullptr);
	return get<double>(mp.m_index);
}

vec3d FENodeDataMap::valueVec3d(const FEMaterialPoint& mp)
{
	assert(mp.m_elem == nullptr);
	return get<vec3d>(mp.m_index);
}

mat3d FENodeDataMap::valueMat3d(const FEMaterialPoint& mp)
{
	assert(mp.m_elem == nullptr);
	return get<mat3d>(mp.m_index);
}

mat3ds FENodeDataMap::valueMat3ds(const FEMaterialPoint& mp)
{
	assert(mp.m_elem == nullptr);
	return get<mat3ds>(mp.m_index);
}

// return the item list associated with this map
FEItemList* FENodeDataMap::GetItemList()
{
	return const_cast<FENodeSet*>(m_nodeSet);
}

void FENodeDataMap::Serialize(DumpStream& ar)
{
	FEDataMap::Serialize(ar);
	if (ar.IsShallow() == false)
	{
		if (ar.IsSaving())
		{
			// We have to cast the const away before serializing
			FENodeSet* ns = const_cast<FENodeSet*>(m_nodeSet);
			ar << ns;
		}
		else
		{
			FENodeSet* ns;
			ar >> ns;
			m_nodeSet = ns;
		}
	}
}
