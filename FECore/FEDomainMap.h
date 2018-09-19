#pragma once
#include "FEDataMap.h"
#include "FEDomain.h"
#include "FEElementSet.h"
#include "FEMaterialPoint.h"

class FEDomainMap : public FEDataMap
{
public:
	//! default constructor
	FEDomainMap(int dataType);

	//! copy constructor
	FEDomainMap(const FEDomainMap& map);

	//! assignment operator
	FEDomainMap& operator = (const FEDomainMap& map);

	//! Create a surface data map for this surface
	bool Create(const FEDomain* ps, double val = 0.0);

	//! Create a surface data map for this surface
	bool Create(FEElementSet* ps, double val = 0.0);

	//! serialization
	void Serialize(DumpStream& ar);

	//! set the name
	void SetName(const std::string& name);

	//! get the name
	const std::string& GetName() const { return m_name; }

	//! get the value at a material point
	double value(const FEMaterialPoint& pt) override;
	vec3d valueVec3d(const FEMaterialPoint& pt) override;

public:
	template <typename T> T value(int nface, int node);
	template <typename T> void setValue(int nface, int node, const T& v);

	void setValue(int n, double v);
	void setValue(int n, const vec2d& v);
	void setValue(int n, const vec3d& v);

	void fillValue(double v);
	void fillValue(const vec2d& v);
	void fillValue(const vec3d& v);

private:
	int	m_maxElemNodes;	// number of nodes for each element
	std::string	m_name;
	const FEDomain*	m_dom;
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
