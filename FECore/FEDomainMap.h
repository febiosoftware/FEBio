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
	mat3ds valueMat3ds(const FEMaterialPoint& pt) override;

	//! Get the element set
	const FEElementSet* GetElementSet() const { return m_elset; }

	//! return max nr of nodes
	int MaxNodes() const { return m_maxElemNodes; }

	//! return storage format
	int	StorageFormat() const { return m_fmt; }

	// return the item list associated with this map
	FEItemList* GetItemList() override;

	// merge with another map
	bool Merge(FEDomainMap& map);

public:
	template <typename T> T value(int nelem, int node)
	{
		return get<T>(nelem*m_maxElemNodes + node);
	}

	template <typename T> void setValue(int nelem, int node, const T& v)
	{
		set<T>(nelem*m_maxElemNodes + node, v);
	}

	void setValue(int n, double v) override;
	void setValue(int n, const vec2d& v) override;
	void setValue(int n, const vec3d& v) override;
	void setValue(int n, const mat3d& v) override;
	void setValue(int n, const mat3ds& v) override;

	void fillValue(double v) override;
	void fillValue(const vec2d& v) override;
	void fillValue(const vec3d& v) override;
	void fillValue(const mat3d& v) override;
	void fillValue(const mat3ds& v) override;

private:
	void Realloc(int newElemSize, int newMaxElemNodes);

private:
	int					m_fmt;				//!< storage format
	int					m_maxElemNodes;		//!< max number of nodes for each element
	FEElementSet*		m_elset;			//!< the element set on which this map is defined

	vector<int>		m_NLT;		//!< node index lookup table for FMT_NODE
	int				m_imin;		//!< min index for lookup for FMT_NODE
};
