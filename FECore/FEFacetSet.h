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
#include "fecore_api.h"
#include "FEItemList.h"
#include "FEElement.h"
#include "FENodeSet.h"
#include <vector>
#include <string>

class FESurface;

//-----------------------------------------------------------------------------
//! This class defines a set of facets. This can be used in the creation of
//! surfaces.
class FECORE_API FEFacetSet : public FEItemList
{
public:
	struct FACET
	{
		// max nr of nodes for each facet
		enum { MAX_NODES = 10};

		// different facet types
		enum FacetType {
			INVALID = 0,
			TRI3  = 3,
			QUAD4 = 4,
			TRI6  = 6,
			TRI7  = 7,
			QUAD8 = 8,
			QUAD9 = 9,
			TRI10 = 10
		};

		int	node[FACET::MAX_NODES];
		int	ntype;	//	3=tri3, 4=quad4, 6=tri6, 7=tri7, 8=quad8, 9=quad9

		void Serialize(DumpStream& ar);

		FACET() { ntype = FACET::INVALID; }
	};

public:
	// Constructor
	FEFacetSet(FEModel* fem);

	// Allocate facets 
	void Create(int n);

	// create from a surface
	void Create(const FESurface& surf);

	// return the size of the facet ste
	int Faces() const;

	// return a facet
	FACET& Face(int i);
	const FACET& Face(int i) const;

	// add a facet set
	void Add(FEFacetSet* pf);

	// extract the node set from the facet set
	FENodeList GetNodeList() const;

public:
	// serialize
	void Serialize(DumpStream& ar);

	static void SaveClass(DumpStream& ar, FEFacetSet* p);
	static FEFacetSet* LoadClass(DumpStream& ar, FEFacetSet* p);

	// TODO: This is a hack used to convert between a surface and a facet set when reading face data records
	void SetSurface(FESurface* surf);
	FESurface* GetSurface();

private:
	std::vector<FACET>	m_Face;	// the list of facets
	FESurface*			m_surface;		// the surface this facet list refers to (can be zero!)
};
