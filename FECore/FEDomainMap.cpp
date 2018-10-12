#include "stdafx.h"
#include "FEDomainMap.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FEDomainMap::FEDomainMap(FEDataType dataType) : FEDataMap(dataType)
{
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
	int NF = ps->Elements();
	FEMesh* mesh = ps->GetMesh();
	m_maxElemNodes = 0;
	for (int i = 0; i<NF; ++i)
	{
		FEElement& el = *mesh->FindElementFromID((*ps)[i]);
		int ne = el.Nodes();
		if (ne > m_maxElemNodes) m_maxElemNodes = ne;
	}
	return resize(NF*m_maxElemNodes*DataSize(), val);
}

//-----------------------------------------------------------------------------
void FEDomainMap::setValue(int n, double v)
{
	int index = n*m_maxElemNodes;
	for (int i = 0; i<m_maxElemNodes; ++i) set<double>(index + i, v);
}

//-----------------------------------------------------------------------------
void FEDomainMap::setValue(int n, const vec2d& v)
{
	int index = n*m_maxElemNodes;
	for (int i = 0; i<m_maxElemNodes; ++i) set<vec2d>(index + i, v);
}

//-----------------------------------------------------------------------------
void FEDomainMap::setValue(int n, const vec3d& v)
{
	int index = n*m_maxElemNodes;
	for (int i = 0; i<m_maxElemNodes; ++i) set<vec3d>(index + i, v);
}

//-----------------------------------------------------------------------------
void FEDomainMap::setValue(int n, const mat3d& v)
{
	int index = n*m_maxElemNodes;
	for (int i = 0; i<m_maxElemNodes; ++i) set<mat3d>(index + i, v);
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
	if (ar.IsSaving())
	{
		ar << m_maxElemNodes;
		ar << m_name;
	}
	else
	{
		ar >> m_maxElemNodes;
		ar >> m_name;
	}
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

	// get shape functions
	double* H = pe->H(pt.m_index);

	double v = 0.0;
	int ne = pe->Nodes();
	for (int i = 0; i < ne; ++i)
	{
		double vi = value<double>(lid, i);
		v += vi*H[i];
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

	// get shape functions
	double* H = pe->H(pt.m_index);

	vec3d v(0,0,0);
	int ne = pe->Nodes();
	for (int i = 0; i < ne; ++i)
	{
		vec3d vi = value<vec3d>(lid, i);
		v += vi*H[i];
	}

	return v;
}
