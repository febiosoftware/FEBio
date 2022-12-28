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
#include "FEElement.h"
#include <vector>
#include <string>
#include "FEItemList.h"
#include "FENodeList.h"

//-----------------------------------------------------------------------------
//! This class defines a set of segments. This can be used in the creation of edges.
class FECORE_API FESegmentSet : public FEItemList
{
public:
	struct SEGMENT
	{
		enum SegmentType {
			INVALID=0,
			LINE2 = 2,
			LINE3 = 3
		};

		int	node[3];
		int	ntype;	//	2=line2

		void Serialize(DumpStream& ar);

		SEGMENT() { ntype = SEGMENT::INVALID; }
	};

public:
	// constructor
	FESegmentSet(FEModel* fem);

	// allocate segments
	void Create(int n);

	// return nr of segments
	int Segments() const { return (int)m_Seg.size(); }

	// return a segment
	SEGMENT& Segment(int i);
	const SEGMENT& Segment(int i) const;

	FENodeList GetNodeList() const;

public:
	// serialization
	void Serialize(DumpStream& ar);

	static void SaveClass(DumpStream& ar, FESegmentSet* p);
	static FESegmentSet* LoadClass(DumpStream& ar, FESegmentSet* p);

private:
	vector<SEGMENT>	m_Seg;	// the actual segment list
};
