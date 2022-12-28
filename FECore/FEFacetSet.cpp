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
#include "FEFacetSet.h"
#include "FEMesh.h"
#include "DumpStream.h"
#include "FESurface.h"

//-----------------------------------------------------------------------------
void FEFacetSet::FACET::Serialize(DumpStream& ar)
{
	ar & node;
	ar & ntype;
}

//-----------------------------------------------------------------------------
FEFacetSet::FEFacetSet(FEModel* fem) : FEItemList(fem)
{
	m_surface = nullptr;
}

//-----------------------------------------------------------------------------
void FEFacetSet::SetSurface(FESurface* surf) { m_surface = surf; }
FESurface* FEFacetSet::GetSurface() { return m_surface; }

//-----------------------------------------------------------------------------
int FEFacetSet::Faces() const { return (int)m_Face.size(); }

//-----------------------------------------------------------------------------
void FEFacetSet::Create(int n)
{
	m_Face.resize(n);
}

//-----------------------------------------------------------------------------
// create from a surface
void FEFacetSet::Create(const FESurface& surf)
{
	int NE = surf.Elements();
	m_Face.resize(NE);
	for (int i = 0; i < NE; ++i)
	{
		const FESurfaceElement& el = surf.Element(i);
		FACET& face = m_Face[i];
		switch (el.Shape())
		{
		case ET_TRI3 : face.ntype = 3; break;
		case ET_QUAD4: face.ntype = 4; break;
		case ET_TRI6 : face.ntype = 6; break;
		case ET_TRI7 : face.ntype = 7; break;
		case ET_QUAD8: face.ntype = 8; break;
		default:
			assert(false);
		}

		for (int j = 0; j < el.Nodes(); ++j)
		{
			face.node[j] = el.m_node[j];
		}
	}
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

void FEFacetSet::SaveClass(DumpStream& ar, FEFacetSet* p) {}
FEFacetSet* FEFacetSet::LoadClass(DumpStream& ar, FEFacetSet* p)
{
	p = new FEFacetSet(&ar.GetFEModel());
	return p;
}
