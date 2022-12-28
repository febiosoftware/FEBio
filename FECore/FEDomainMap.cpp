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
#include "FEDomainMap.h"
#include "FEMesh.h"
#include "DumpStream.h"
#include "mathalg.h"

//-----------------------------------------------------------------------------
FEDomainMap::FEDomainMap() : FEDataMap(FE_DOMAIN_MAP)
{
	m_maxElemNodes = 0;
}

//-----------------------------------------------------------------------------
FEDomainMap::FEDomainMap(FEDataType dataType, Storage_Fmt format) : FEDataMap(FE_DOMAIN_MAP, dataType)
{
	m_fmt = format;
	m_maxElemNodes = 0;
}

//-----------------------------------------------------------------------------
FEDomainMap::FEDomainMap(const FEDomainMap& map) : FEDataMap(map)
{
	m_name = map.m_name;
	m_maxElemNodes = map.m_maxElemNodes;
}

//-----------------------------------------------------------------------------
FEDomainMap& FEDomainMap::operator = (const FEDomainMap& map)
{
	FEDataArray::operator=(map);
	m_name = map.m_name;
	m_maxElemNodes = map.m_maxElemNodes;
	return *this;
}

//-----------------------------------------------------------------------------
// return the item list associated with this map
FEItemList* FEDomainMap::GetItemList()
{
	return m_elset;
}

//-----------------------------------------------------------------------------
bool FEDomainMap::Create(FEElementSet* ps, double val)
{
	m_elset = ps;
	int NE = ps->Elements();
	FEMesh* mesh = ps->GetMesh();
	m_maxElemNodes = 0;

	if (m_fmt == FMT_MULT)
	{
		for (int i = 0; i < NE; ++i)
		{
			FEElement& el = *mesh->FindElementFromID((*ps)[i]);
			int ne = el.Nodes();
			if (ne > m_maxElemNodes) m_maxElemNodes = ne;
		}
		return resize(NE*m_maxElemNodes, val);
	}
	else if (m_fmt == FMT_ITEM)
	{
		m_maxElemNodes = 1;
		return resize(NE, val);
	}
	else if (m_fmt == FMT_NODE)
	{
		FENodeList nodeList = ps->GetNodeList();
		int nodes = nodeList.Size();

		// find min, max index
		m_imin = mesh->Nodes();
		int imax = -1;
		for (int i = 0; i < nodes; ++i)
		{
			int nid = nodeList[i];
			if (nid < m_imin) m_imin = nid;
			if (nid > imax) imax = nid;
		}

		// build lookup table
		m_NLT.resize(imax - m_imin + 1, -1);
		for (int i = 0; i < nodes; ++i)
		{
			int nid = nodeList[i];
			m_NLT[nid - m_imin] = i;
		}
		
		return resize(nodes, val);
	}
	else if (m_fmt == FMT_MATPOINTS)
	{
		for (int i = 0; i < NE; ++i)
		{
			FEElement& el = *mesh->FindElementFromID((*ps)[i]);
			int ne = el.GaussPoints();
			if (ne > m_maxElemNodes) m_maxElemNodes = ne;
		}
		return resize(NE*m_maxElemNodes, val);
	}
	else return false;
}

//-----------------------------------------------------------------------------
void FEDomainMap::setValue(int n, double v)
{
	if (m_fmt == FMT_MULT)
	{
		int index = n*m_maxElemNodes;
		for (int i = 0; i < m_maxElemNodes; ++i) set<double>(index + i, v);
	}
	else if (m_fmt == FMT_ITEM)
	{
		set<double>(n, v);
	}
	else if (m_fmt == FMT_NODE)
	{
		set<double>(n, v);
	}
}

//-----------------------------------------------------------------------------
void FEDomainMap::setValue(int n, const vec2d& v)
{
	if (m_fmt == FMT_MULT)
	{
		int index = n*m_maxElemNodes;
		for (int i = 0; i < m_maxElemNodes; ++i) set<vec2d>(index + i, v);
	}
	else if (m_fmt == FMT_ITEM)
	{
		set<vec2d>(n, v);
	}
	else if (m_fmt == FMT_NODE)
	{
		set<vec2d>(n, v);
	}
}

//-----------------------------------------------------------------------------
void FEDomainMap::setValue(int n, const vec3d& v)
{
	if ((m_fmt == FMT_MULT) || (m_fmt == FMT_MATPOINTS))
	{
		int index = n*m_maxElemNodes;
		for (int i = 0; i < m_maxElemNodes; ++i) set<vec3d>(index + i, v);
	}
	else if (m_fmt == FMT_ITEM)
	{
		set<vec3d>(n, v);
	}
	else if (m_fmt == FMT_NODE)
	{
		set<vec3d>(n, v);
	}
}

//-----------------------------------------------------------------------------
void FEDomainMap::setValue(int n, const mat3d& v)
{
	if ((m_fmt == FMT_MULT) || (m_fmt == FMT_MATPOINTS))
	{
		int index = n*m_maxElemNodes;
		for (int i = 0; i < m_maxElemNodes; ++i) set<mat3d>(index + i, v);
	}
	else if (m_fmt == FMT_ITEM)
	{
		set<mat3d>(n, v);
	}
	else if (m_fmt == FMT_NODE)
	{
		set<mat3d>(n, v);
	}
}

//-----------------------------------------------------------------------------
void FEDomainMap::setValue(int n, const mat3ds& v)
{
	if ((m_fmt == FMT_MULT) || (m_fmt == FMT_MATPOINTS))
	{
		int index = n * m_maxElemNodes;
		for (int i = 0; i < m_maxElemNodes; ++i) set<mat3ds>(index + i, v);
	}
	else if (m_fmt == FMT_ITEM)
	{
		set<mat3ds>(n, v);
	}
	else if (m_fmt == FMT_NODE)
	{
		set<mat3ds>(n, v);
	}
}

//-----------------------------------------------------------------------------
void FEDomainMap::fillValue(double v)
{
	set<double>(v);
}

//-----------------------------------------------------------------------------
void FEDomainMap::fillValue(const vec2d& v)
{
	set<vec2d>(v);
}

//-----------------------------------------------------------------------------
void FEDomainMap::fillValue(const vec3d& v)
{
	set<vec3d>(v);
}

//-----------------------------------------------------------------------------
void FEDomainMap::fillValue(const mat3d& v)
{
	set<mat3d>(v);
}

//-----------------------------------------------------------------------------
void FEDomainMap::fillValue(const mat3ds& v)
{
	set<mat3ds>(v);
}

//-----------------------------------------------------------------------------
void FEDomainMap::Serialize(DumpStream& ar)
{
	FEDataArray::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_maxElemNodes;
	ar & m_name;
	ar & m_elset;
	ar & m_fmt;
}

//-----------------------------------------------------------------------------
//! get the value at a material point
double FEDomainMap::value(const FEMaterialPoint& pt)
{
	// get the element this material point is in
	FEElement* pe = pt.m_elem;
	assert(pe);

	// see if this element belongs to the element set
	assert(m_elset);
	int lid = m_elset->GetLocalIndex(*pe);
	assert((lid >= 0));

	double v = 0.0;
	if (m_fmt == FMT_MULT)
	{
		// get shape functions
		double* H = pe->H(pt.m_index);

		int ne = pe->Nodes();
		for (int i = 0; i < ne; ++i)
		{
			double vi = value<double>(lid, i);
			v += vi*H[i];
		}
	}
	else if (m_fmt == FMT_ITEM)
	{
		v = get<double>(lid);
	}
	else if (m_fmt == FMT_MATPOINTS)
	{
		v = value<double>(lid, pt.m_index);
	}
	else { assert(false); }

	return v;
}

//-----------------------------------------------------------------------------
//! get the value at a material point
vec3d FEDomainMap::valueVec3d(const FEMaterialPoint& pt)
{
	// get the element this material point is in
	FEElement* pe = pt.m_elem;
	assert(pe);

	// see if this element belongs to the element set
	assert(m_elset);
	int lid = m_elset->GetLocalIndex(*pe);
	assert((lid >= 0));

	vec3d v(0, 0, 0);
	if (m_fmt == FMT_MULT)
	{
		// get shape functions
		double* H = pe->H(pt.m_index);
		int ne = pe->Nodes();
		for (int i = 0; i < ne; ++i)
		{
			vec3d vi = value<vec3d>(lid, i);
			v += vi*H[i];
		}
	}
	else if (m_fmt == FMT_ITEM)
	{
		return get<vec3d>(lid);
	}

	return v;
}

//-----------------------------------------------------------------------------
//! get the value at a material point
mat3d FEDomainMap::valueMat3d(const FEMaterialPoint& pt)
{
	assert(DataType() == FEDataType::FE_MAT3D);
	// get the element this material point is in
	FEElement* pe = pt.m_elem;
	assert(pe);

	// see if this element belongs to the element set
	assert(m_elset);
	int lid = m_elset->GetLocalIndex(*pe);
	assert((lid >= 0));

	mat3d Q;
	if (m_fmt == FMT_ITEM)
	{
		Q = get<mat3d>(lid);
	}
	else if (m_fmt == FMT_MATPOINTS)
	{
		Q = value<mat3d>(lid, pt.m_index);
	}
	else { assert(false); }

	return Q;
}

//-----------------------------------------------------------------------------
//! get the value at a material point
mat3ds FEDomainMap::valueMat3ds(const FEMaterialPoint& pt)
{
	assert(DataType() == FEDataType::FE_MAT3DS);
	// get the element this material point is in
	FEElement* pe = pt.m_elem;
	assert(pe);

	// see if this element belongs to the element set
	assert(m_elset);
	int lid = m_elset->GetLocalIndex(*pe);
	assert((lid >= 0));

	mat3ds Q;
	if (m_fmt == FMT_ITEM)
	{
		Q = get<mat3ds>(lid);
	}
	else if (m_fmt == FMT_NODE)
	{
		// calculate shape function values
		double* w = pe->H(pt.m_index);
		mat3ds Qi[FEElement::MAX_NODES];
		int ne = pe->Nodes();
		for (int i = 0; i < ne; ++i)
		{
			int nid = pe->m_node[i];

			int lid = m_NLT[nid] - m_imin; 
			assert((lid >= 0) && (lid < DataCount()));

			Qi[i] = get<mat3ds>(lid);
		}

		// weighted average
		Q = weightedAverageStructureTensor(Qi, w, ne);
	}
	else if (m_fmt == FMT_MATPOINTS)
	{
		return value<mat3ds>(lid, pt.m_index);
	}

	return Q;
}

//-----------------------------------------------------------------------------
// merge with another map
bool FEDomainMap::Merge(FEDomainMap& map)
{
	// make sure the type and format are the same
	if (map.DataType() != DataType()) return false;
	if (map.StorageFormat() != map.StorageFormat()) return false;

	// get the two element sets
	FEElementSet* set1 = m_elset;
	const FEElementSet* set2 = map.GetElementSet();
	assert(set1->GetMesh() == set2->GetMesh());

	// create a new element set
	// TODO: should we add it to the mesh? I think we probably have to for remeshing
	FEElementSet* elset = new FEElementSet(set1->GetMesh());
	elset->Add(*set1);
	elset->Add(*set2);

	// reallocate the data array
	int oldElems = set1->Elements();
	int newElems = elset->Elements();

	if (StorageFormat() == FMT_MULT)
	{
		assert(DataType() == FE_DOUBLE);

		// get the max element nodes
		int maxElemNodes = m_maxElemNodes;
		if (map.m_maxElemNodes > maxElemNodes) maxElemNodes = map.m_maxElemNodes;
		Realloc(newElems, maxElemNodes);

		// set the new values of the map
		for (int i = 0; i < set2->Elements(); ++i)
		{
			int n = set2->Element(i).Nodes();
			for (int j = 0; j < n; ++j)
			{
				double v = map.value<double>(i, j);
				setValue<double>(oldElems + i, j, v);
			}
		}
	}
	else if (StorageFormat() == FMT_MATPOINTS)
	{
		assert(DataType() == FE_DOUBLE);

		// get the max element nodes
		int maxElemNodes = m_maxElemNodes;
		if (map.m_maxElemNodes > maxElemNodes) maxElemNodes = map.m_maxElemNodes;
		Realloc(newElems, maxElemNodes);

		// set the new values of the map
		for (int i = 0; i < set2->Elements(); ++i)
		{
			int n = set2->Element(i).GaussPoints();
			for (int j = 0; j < n; ++j)
			{
				double v = map.value<double>(i, j);
				setValue<double>(oldElems + i, j, v);
			}
		}
	}
	else if (StorageFormat() == FMT_ITEM)
	{
		realloc(newElems);

		switch (DataType())
		{
		case FE_DOUBLE:
		{
			// set the new values of the map
			for (int i = 0; i < set2->Elements(); ++i)
			{
				double v = map.get<double>(i);
				set<double>(oldElems + i, v);
			}
		}
		break;
		case FE_MAT3D:
		{
			// set the new values of the map
			for (int i = 0; i < set2->Elements(); ++i)
			{
				mat3d v = map.get<mat3d>(i);
				set<mat3d>(oldElems + i, v);
			}
		}
		break;
		}
	}

	m_elset = elset;

	return true;
}

//-----------------------------------------------------------------------------
void FEDomainMap::Realloc(int newElemSize, int newMaxElemNodes)
{
	int oldMaxElemNodes = m_maxElemNodes;
	assert(newMaxElemNodes >= oldMaxElemNodes);
	if (newMaxElemNodes == oldMaxElemNodes) realloc(newElemSize*newMaxElemNodes);
	else
	{
		FEDomainMap oldData(*this);
		m_maxElemNodes = newMaxElemNodes;
		realloc(newElemSize*newMaxElemNodes);
		FEElementSet* set = m_elset;
		for (int i = 0; i < set->Elements(); ++i)
		{
			int n = set->Element(i).Nodes();
			if (StorageFormat() == FMT_MATPOINTS) n = set->Element(i).GaussPoints();

			for (int j = 0; j < n; ++j)
			{
				double v = oldData.get<double>(i*oldMaxElemNodes + j);
				setValue<double>(i, j, v);
			}
		}
	}
}
