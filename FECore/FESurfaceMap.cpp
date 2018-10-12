#include "stdafx.h"
#include "FESurfaceMap.h"
#include "FESurface.h"
#include "DumpStream.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap(FEDataType dataType) : FEDataMap(dataType)
{
	m_maxFaceNodes = 0;
}

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap(const FESurfaceMap& map) : FEDataMap(map), m_name(map.m_name)
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
	m_dom = ps;
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
	m_dom = 0;
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
void FESurfaceMap::setValue(int n, const mat3d& v)
{
	int index = n*m_maxFaceNodes;
	for (int i = 0; i<m_maxFaceNodes; ++i) set<mat3d>(index + i, v);
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
void FESurfaceMap::fillValue(const mat3d& v)
{
	set<mat3d>(v);
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

//-----------------------------------------------------------------------------
double FESurfaceMap::value(const FEMaterialPoint& pt)
{
	// get the element this material point is in
	FESurfaceElement* pe = dynamic_cast<FESurfaceElement*>(pt.m_elem);
	assert(pe);

	// make sure this element belongs to this domain
	// TODO: Can't check this if map was created through FEFacetSet
//	assert(pe->GetDomain() == m_dom);

	// get its local ID
	int lid = pe->GetLocalID();

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
vec3d FESurfaceMap::valueVec3d(const FEMaterialPoint& pt)
{
	// get the element this material point is in
	FESurfaceElement* pe = dynamic_cast<FESurfaceElement*>(pt.m_elem);
	assert(pe);

	// make sure this element belongs to this domain
	// TODO: Can't check this if map was created through FEFacetSet
	//	assert(pe->GetDomain() == m_dom);

	// get its local ID
	int lid = pe->GetLocalID();

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

