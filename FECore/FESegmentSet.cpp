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



#include "stdafx.h"
#include "FESegmentSet.h"
#include "DumpStream.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
void FESegmentSet::SEGMENT::Serialize(DumpStream& ar)
{
	ar & node;
	ar & ntype;
}

//-----------------------------------------------------------------------------
FESegmentSet::FESegmentSet(FEModel* fem) : FEItemList(fem)
{
}

//-----------------------------------------------------------------------------
void FESegmentSet::Create(int n)
{
	m_Seg.resize(n);
}

//-----------------------------------------------------------------------------
FESegmentSet::SEGMENT& FESegmentSet::Segment(int i)
{
	return m_Seg[i];
}

//-----------------------------------------------------------------------------
const FESegmentSet::SEGMENT& FESegmentSet::Segment(int i) const
{
	return m_Seg[i];
}

//-----------------------------------------------------------------------------
void FESegmentSet::Serialize(DumpStream& ar)
{
	FEItemList::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_Seg;
}

void FESegmentSet::SaveClass(DumpStream& ar, FESegmentSet* p)
{

}

FESegmentSet* FESegmentSet::LoadClass(DumpStream& ar, FESegmentSet* p)
{
	p = new FESegmentSet(&ar.GetFEModel());
	return p;
}

//-----------------------------------------------------------------------------
FENodeList FESegmentSet::GetNodeList() const
{
	FEMesh* mesh = GetMesh();
	FENodeList set(mesh);
	vector<int> tag(mesh->Nodes(), 0);
	for (int i = 0; i < Segments(); ++i)
	{
		const SEGMENT& el = m_Seg[i];
		int ne = el.ntype;
		for (int j = 0; j < ne; ++j)
		{
			if (tag[el.node[j]] == 0)
			{
				set.Add(el.node[j]);
				tag[el.node[j]] = 1;
			}
		}
	}
	return set;
}
