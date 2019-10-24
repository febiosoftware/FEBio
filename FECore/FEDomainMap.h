/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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



#pragma once
#include "FEDataMap.h"
#include "FEDomain.h"
#include "FEElementSet.h"
#include "FEMaterialPoint.h"

class FECORE_API FEDomainMap : public FEDataMap
{
public:
	//! default constructor
	FEDomainMap();
	FEDomainMap(FEDataType dataType, Storage_Fmt format = FMT_MULT);

	//! copy constructor
	FEDomainMap(const FEDomainMap& map);

	//! assignment operator
	FEDomainMap& operator = (const FEDomainMap& map);

	//! Create a surface data map for this surface
	bool Create(FEElementSet* ps, double val = 0.0);

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! get the value at a material point
	double value(const FEMaterialPoint& pt) override;
	vec3d valueVec3d(const FEMaterialPoint& pt) override;
	mat3d valueMat3d(const FEMaterialPoint& pt) override;

	//! Get the element set
	const FEElementSet* GetElementSet() const { return m_elset; }

	//! return max nr of nodes
	int MaxNodes() const { return m_maxElemNodes; }

	//! return storage format
	int	StorageFormat() const { return m_fmt; }

	// return the item list associated with this map
	FEItemList* GetItemList() override;

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
	int					m_fmt;				//!< storage format
	int					m_maxElemNodes;		//!< max number of nodes for each element
	FEElementSet*		m_elset;			//!< the element set on which this map is defined
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
