#include "stdafx.h"
#include "FEDomainMap.h"
#include "FEMesh.h"
#include "DumpStream.h"

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
		return resize(NE, val);
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
}

//-----------------------------------------------------------------------------
void FEDomainMap::setValue(int n, const vec3d& v)
{
	if (m_fmt == FMT_MULT)
	{
		int index = n*m_maxElemNodes;
		for (int i = 0; i < m_maxElemNodes; ++i) set<vec3d>(index + i, v);
	}
	else if (m_fmt == FMT_ITEM)
	{
		set<vec3d>(n, v);
	}
}

//-----------------------------------------------------------------------------
void FEDomainMap::setValue(int n, const mat3d& v)
{
	if (m_fmt == FMT_MULT)
	{
		int index = n*m_maxElemNodes;
		for (int i = 0; i < m_maxElemNodes; ++i) set<mat3d>(index + i, v);
	}
	else if (m_fmt == FMT_ITEM)
	{
		set<mat3d>(n, v);
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
void FEDomainMap::Serialize(DumpStream& ar)
{
	FEDataArray::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_maxElemNodes;
	ar & m_name;
	ar & m_elset;
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

	return Q;
}
