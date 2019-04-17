/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEFacetSet.h"
#include "FEMesh.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
void FEFacetSet::FACET::Serialize(DumpStream& ar)
{
	ar & node;
	ar & ntype;
}

//-----------------------------------------------------------------------------
FEFacetSet::FEFacetSet(FEModel* fem) : FEItemList(fem)
{
}

//-----------------------------------------------------------------------------
void FEFacetSet::Create(int n)
{
	m_Face.resize(n);
}

//-----------------------------------------------------------------------------
FEFacetSet::FACET& FEFacetSet::Face(int i)
{
	return m_Face[i];
}

//-----------------------------------------------------------------------------
const FEFacetSet::FACET& FEFacetSet::Face(int i) const
{
	return m_Face[i];
}

//-----------------------------------------------------------------------------
void FEFacetSet::Add(FEFacetSet* pf)
{
	m_Face.insert(m_Face.end(), pf->m_Face.begin(), pf->m_Face.end());
}

//-----------------------------------------------------------------------------
FENodeList FEFacetSet::GetNodeList() const
{
	FEMesh* mesh = GetMesh();
	FENodeList set(mesh);
	vector<int> tag(mesh->Nodes(), 0);
	for (int i = 0; i<Faces(); ++i)
	{
		const FACET& el = m_Face[i];
		int ne = el.ntype;
		for (int j = 0; j<ne; ++j)
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

//-----------------------------------------------------------------------------
void FEFacetSet::Serialize(DumpStream& ar)
{
	FEItemList::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_Face;
}
