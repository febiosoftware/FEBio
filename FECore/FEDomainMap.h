#pragma once
#include "FEDataMap.h"
#include "FEDomain.h"
#include "FEElementSet.h"
#include "FEMaterialPoint.h"

class FEDomainMap : public FEDataMap
{
public:
	//! default constructor
	FEDomainMap(FEDataType dataType, Storage_Fmt format = FMT_MULT);

	//! copy constructor
	FEDomainMap(const FEDomainMap& map);

	//! assignment operator
	FEDomainMap& operator = (const FEDomainMap& map);

	//! Create a surface data map for this surface
	bool Create(FEElementSet* ps, double val = 0.0);

	//! serialization
	void Serialize(DumpStream& ar);

	//! get the value at a material point
	double value(const FEMaterialPoint& pt) override;
	vec3d valueVec3d(const FEMaterialPoint& pt) override;
	mat3d valueMat3d(const FEMaterialPoint& pt) override;

	//! Get the element set
	const FEElementSet* GetElementSet() const { return m_elset; }

	//! return max nr of nodes
	int MaxNodes() const { return m_maxElemNodes; }

	//! return storage format
	Storage_Fmt	StorageFormat() const { return m_fmt; }

public:
	template <typename T> T value(int nface, int node);
	template <typename T> void setValue(int nface, int node, const T& v);

	void setValue(int n, double v) override;
	void setValue(int n, const vec2d& v) override;
	void setValue(int n, const vec3d& v) override;
	void setValue(int n, const mat3d& v) override;

	void fillValue(double v) override;
	void fillValue(const vec2d& v) override;
	void fillValue(const vec3d& v) override;
	void fillValue(const mat3d& v) override;

private:
	Storage_Fmt				m_fmt;				//!< storage format
	int						m_maxElemNodes;		//!< max number of nodes for each element
	const FEElementSet*		m_elset;			//!< the element set on which this map is defined
};

template <> inline double FEDomainMap::value(int nelem, int node)
{
	return get<double>(nelem*m_maxElemNodes + node);
}

template <> inline vec2d FEDomainMap::value(int nelem, int node)
{
	return get<vec2d>(nelem*m_maxElemNodes + node);
}

template <> inline vec3d FEDomainMap::value(int nelem, int node)
{
	return get<vec3d>(nelem*m_maxElemNodes + node);
}

template <> inline void FEDomainMap::setValue(int nelem, int node, const double& v)
{
	set<double>(nelem*m_maxElemNodes + node, v);
}

template <> inline void FEDomainMap::setValue(int nelem, int node, const vec2d& v)
{
	set<vec2d>(nelem*m_maxElemNodes + node, v);
}

template <> inline void FEDomainMap::setValue(int nelem, int node, const vec3d& v)
{
	set<vec3d>(nelem*m_maxElemNodes + node, v);
}

template <> inline void FEDomainMap::setValue(int nelem, int node, const mat3d& v)
{
	set<mat3d>(nelem*m_maxElemNodes + node, v);
}
