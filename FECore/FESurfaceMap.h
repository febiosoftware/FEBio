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
#include <vector>
#include <string>
#include <assert.h>
#include "FEDataMap.h"

//-----------------------------------------------------------------------------
class FESurface;
class FEFacetSet;
class DumpStream;

//-----------------------------------------------------------------------------
typedef int FEFacetIndex;


//-----------------------------------------------------------------------------
// TODO: Perhaps I should rename this FESurfaceData
//       and then define FESurfaceMap as a tool for evaluating data across a surface (i.e. via shape functions)
class FECORE_API FESurfaceMap : public FEDataMap
{
public:
	//! default constructor
	FESurfaceMap();
	FESurfaceMap(FEDataType dataType);

	//! copy constructor
	FESurfaceMap(const FESurfaceMap& map);

	//! assignment operator
	FESurfaceMap& operator = (const FESurfaceMap& map);

	//! Create a surface data map for this surface
	bool Create(const FEFacetSet* ps, double val = 0.0, Storage_Fmt fmt = FMT_MULT);

	//! serialization
	void Serialize(DumpStream& ar) override;

	const FEFacetSet* GetFacetSet() const { return m_surf; }

	int MaxNodes() const { return m_maxFaceNodes; }

	// return the item list associated with this map
	FEItemList* GetItemList() override;

	int StorageFormat() const;

public: // from FEDataMap
	double value(const FEMaterialPoint& pt) override;
	vec3d valueVec3d(const FEMaterialPoint& mp) override;
	mat3d valueMat3d(const FEMaterialPoint& mp) override;
	mat3ds valueMat3ds(const FEMaterialPoint& mp) override;

public:
	template <typename T> T value(int nface, int node)
	{
		return get<T>(nface*m_maxFaceNodes + node);
	}

	template <typename T> void setValue(int nface, int node, const T& v)
	{
		set<T>(nface*m_maxFaceNodes + node, v);
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
	const FEFacetSet*	m_surf;		// the surface for which this data set is defined
	int					m_format;	// the storage format
	int	m_maxFaceNodes;				// number of nodes for each face
};
