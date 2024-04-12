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
#include "FEMeshPartition.h"
#include "FENodeList.h"
#include "FESegmentSet.h"
#include <vector>

//-----------------------------------------------------------------------------
class FECORE_API FELineMaterialPoint : public FEMaterialPoint
{
public:
	// return the surface element
	FELineElement* LineElement() { return (FELineElement*)m_elem; }

	void Serialize(DumpStream& ar) override
	{
		FEMaterialPoint::Serialize(ar);
	}
};

//-----------------------------------------------------------------------------
// This class represents an edge of a domain.
class FECORE_API FEEdge : public FEMeshPartition
{
public:
	FECORE_SUPER_CLASS(FEEDGE_ID)
	FECORE_BASE_CLASS(FEEdge)

public:
	//! constructor
	FEEdge(FEModel* fem);

	//! destructor
	virtual ~FEEdge();

	//! initialize edge data structure
	virtual bool Init() override;

	//! creates edge
	void Create(int nelems, int elemType = -1);

	//! create from edge set
	virtual bool Create(FESegmentSet& es);

	//! extract node set
	FENodeList GetNodeList();

public:

	//! return number of edge elements
	int Elements() const override { return (int)m_Elem.size(); }

	//! return an element of the edge
	FELineElement& Element(int i) { return m_Elem[i]; }

	//! returns reference to element
	FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

	// Create material point data for this surface
	virtual FEMaterialPoint* CreateMaterialPoint();

	void CreateMaterialPointData();

public:
	// Get current coordinates
	void GetNodalCoordinates(FELineElement& el, vec3d* rt);

	// Get reference coordinates
	void GetReferenceNodalCoordinates(FELineElement& el, vec3d* rt);

protected:
	bool Create(FESegmentSet& es, int elemType);

protected:
	std::vector<FELineElement>	m_Elem;
};
