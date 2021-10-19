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
#include "FEFaceList.h"
#include "FEMesh.h"
#include "FEDomain.h"
#include "FEElemElemList.h"
#include "FEElementList.h"
#include "FEEdgeList.h"

bool FEFaceList::FACE::IsEqual(int* n) const
{
	if (ntype == 3)
	{
		if ((node[0] != n[0]) && (node[0] != n[1]) && (node[0] != n[2])) return false;
		if ((node[1] != n[0]) && (node[1] != n[1]) && (node[1] != n[2])) return false;
		if ((node[2] != n[0]) && (node[2] != n[1]) && (node[2] != n[2])) return false;
	}
	else if (ntype == 4)
	{
		if ((node[0] != n[0]) && (node[0] != n[1]) && (node[0] != n[2]) && (node[0] != n[3])) return false;
		if ((node[1] != n[0]) && (node[1] != n[1]) && (node[1] != n[2]) && (node[1] != n[3])) return false;
		if ((node[2] != n[0]) && (node[2] != n[1]) && (node[2] != n[2]) && (node[2] != n[3])) return false;
		if ((node[3] != n[0]) && (node[3] != n[1]) && (node[3] != n[2]) && (node[3] != n[3])) return false;
	}
	else 
	{
		assert(false);
		return false;
	}

	return true;
}

bool FEFaceList::FACE::HasEdge(int a, int b) const
{
	const int* n = node;
	if (ntype == 3)
	{
		if (((n[0] == a) && (n[1] == b)) || ((n[0] == b) && (n[1] == a))) return true;
		if (((n[1] == a) && (n[2] == b)) || ((n[1] == b) && (n[2] == a))) return true;
		if (((n[2] == a) && (n[0] == b)) || ((n[2] == b) && (n[0] == a))) return true;
	}
	else if (ntype == 4)
	{
		if (((n[0] == a) && (n[1] == b)) || ((n[0] == b) && (n[1] == a))) return true;
		if (((n[1] == a) && (n[2] == b)) || ((n[1] == b) && (n[2] == a))) return true;
		if (((n[2] == a) && (n[3] == b)) || ((n[2] == b) && (n[3] == a))) return true;
		if (((n[3] == a) && (n[0] == b)) || ((n[3] == b) && (n[0] == a))) return true;
	}
	else
	{
		assert(false);
	}

	return false;
}

FEFaceList::FEFaceList() : m_mesh(nullptr)
{

}

FEFaceList::FEFaceList(const FEFaceList& faceList)
{
	m_mesh = faceList.m_mesh;
	m_faceList = faceList.m_faceList;
}

int FEFaceList::Faces() const
{
	return (int) m_faceList.size();
}

const FEFaceList::FACE& FEFaceList::operator [] (int i) const
{ 
	return m_faceList[i];
}

const FEFaceList::FACE& FEFaceList::Face(int i) const
{
	return m_faceList[i];
}

FEMesh* FEFaceList::GetMesh()
{
	return m_mesh;
}

int FEElementFaceList::Faces(int elem) const
{
	return (int)m_EFL[elem].size();
}

const std::vector<int>& FEElementFaceList::FaceList(int elem) const
{
	return m_EFL[elem];
}

// Extract the surface only
FEFaceList FEFaceList::GetSurface() const
{
	FEFaceList surface;
	surface.m_mesh = m_mesh;

	int faces = Faces();
	for (int i = 0; i < faces; ++i)
	{
		const FEFaceList::FACE& f = m_faceList[i];
		if (f.nsurf == 1) surface.m_faceList.push_back(f);
	}

	return surface;
}

bool FEFaceList::Create(FEMesh& mesh, FEElemElemList& EEL)
{
	m_mesh = &mesh;

	// get the number of elements in this mesh
	int NE = mesh.Elements();

	// count the number of facets we have to create
	int NF = 0;
	FEElementList EL(mesh);
	FEElementList::iterator it = EL.begin();
	for (int i = 0; i<NE; ++i, ++it)
	{
		FEElement& el = *it;
		int nf = el.Faces();
		for (int j = 0; j<nf; ++j)
		{
			FEElement* pen = EEL.Neighbor(i, j);
			if (pen == 0) ++NF;
			if ((pen != 0) && (el.GetID() < pen->GetID())) ++NF;
		}
	}

	// create the facet list
	m_faceList.resize(NF);

	// build the facets
	int face[FEElement::MAX_NODES];
	NF = 0;
	it = EL.begin();
	for (int i = 0; i<NE; ++i, ++it)
	{
		FEElement& el = *it;
		int nf = el.Faces();
		for (int j = 0; j<nf; ++j)
		{
			FEElement* pen = EEL.Neighbor(i, j);
			if ((pen == 0) || ((pen != 0) && (el.GetID() < pen->GetID())))
			{
				FACE& se = m_faceList[NF++];
				el.GetFace(j, face);

				switch (el.Shape())
				{
				case ET_HEX8:
					se.ntype = 4;
					break;
				case ET_TET4:
				case ET_TET5:
					se.ntype = 3;
					break;
				default:
					assert(false);
				}

				int nn = se.ntype;
				for (int k = 0; k<nn; ++k)
				{
					se.node[k] = face[k];
				}

				// The facet is a surface facet if the element neighbor is null
				se.nsurf = (pen == 0 ? 1 : 0);
			}
		}
	}

	return true;
}

// build the neighbor list
void FEFaceList::BuildNeighbors()
{
	// build a node-face list
	FENodeFaceList NFL;
	NFL.Create(*this);

	for (int i = 0; i < Faces(); ++i)
	{
		FACE& f = m_faceList[i];
		f.nbr[0] = f.nbr[1] = f.nbr[2] = f.nbr[3] = -1;

		if (f.nsurf == 1)
		{
			int fn = f.ntype;
			for (int j = 0; j < fn; ++j)
			{
				int n0 = f.node[j];
				int n1 = f.node[(j + 1) % fn];

				int nval = NFL.Faces(n0);
				std::vector<int> fl = NFL.FaceList(n0);
				for (int k = 0; k < nval; ++k)
				{
					if (fl[k] != i)
					{
						FACE& fk = m_faceList[fl[k]];
						if ((fk.nsurf == 1) && fk.HasEdge(n0, n1))
						{
							f.nbr[j] = fl[k];
							break;
						}
					}
				}
				assert(f.nbr[j] != -1);
			}
		}
	}
}

//=============================================================================
FENodeFaceList::FENodeFaceList()
{

}

int FENodeFaceList::Faces(int node) const
{
	return (int)m_NFL[node].size();
}

const std::vector<int>& FENodeFaceList::FaceList(int node) const
{
	return m_NFL[node];
}

bool FENodeFaceList::Create(FEFaceList& FL)
{
	FEMesh& mesh = *FL.GetMesh();

	int NN = mesh.Nodes();
	m_NFL.resize(NN);

	int NF = FL.Faces();
	for (int i = 0; i < NF; ++i)
	{
		const FEFaceList::FACE& face = FL[i];
		for (int j = 0; j < face.ntype; ++j)
		{
			m_NFL[face.node[j]].push_back(i);
		}
	}

	return true;
}

//=============================================================================

FEElementFaceList::FEElementFaceList()
{

}

//NOTE: only works for hex elements
bool FEElementFaceList::Create(FEElementList& elemList, FEFaceList& faceList)
{
	FEMesh& mesh = *faceList.GetMesh();

	const int FTET[4][3] = {
		{ 0, 1, 3},
		{ 1, 2, 3},
		{ 2, 0, 3},
		{ 2, 1, 0}
	};

	const int FHEX[6][4] = { 
		{ 0, 1, 5, 4 },
		{ 1, 2, 6, 5 },
		{ 2, 3, 7, 6 },
		{ 3, 0, 4, 7 },
		{ 3, 2, 1, 0 },
		{ 4, 5, 6, 7 }};

	// build a node face table for FT to facilitate searching
	int NN = mesh.Nodes();
	vector<vector<int> > NFT; NFT.resize(NN);
	for (int i = 0; i<faceList.Faces(); ++i)
	{
		const FEFaceList::FACE& f = faceList.Face(i);
		NFT[f.node[0]].push_back(i);
		NFT[f.node[1]].push_back(i);
		NFT[f.node[2]].push_back(i);
		if (f.ntype == 4) NFT[f.node[3]].push_back(i);
	}

	int NE = mesh.Elements();
	m_EFL.resize(NE);
	int i = 0;
	for (FEElementList::iterator it = elemList.begin(); it != elemList.end(); ++it, ++i)
	{
		const FEElement& el = *it;
		vector<int>& EFLi = m_EFL[i];
		if ((el.Shape() == ET_TET4) || (el.Shape() == ET_TET5))
		{
			EFLi.resize(4);
			for (int j = 0; j < 4; ++j)
			{
				int fj[3] = { el.m_node[FTET[j][0]], el.m_node[FTET[j][1]], el.m_node[FTET[j][2]] };
				EFLi[j] = -1;
				vector<int>& nfi = NFT[fj[0]];
				for (int k = 0; k<(int)nfi.size(); ++k)
				{
					const FEFaceList::FACE& fk = faceList.Face(nfi[k]);
					if (fk.IsEqual(fj))
					{
						EFLi[j] = nfi[k];
						break;
					}
				}
			}
		}
		else if (el.Shape() == FE_Element_Shape::ET_HEX8)
		{
			EFLi.resize(6);
			for (int j = 0; j < 6; ++j)
			{
				int fj[4] = { el.m_node[FHEX[j][0]], el.m_node[FHEX[j][1]], el.m_node[FHEX[j][2]], el.m_node[FHEX[j][3]] };
				EFLi[j] = -1;
				vector<int>& nfi = NFT[fj[0]];
				for (int k = 0; k<(int)nfi.size(); ++k)
				{
					const FEFaceList::FACE& fk = faceList.Face(nfi[k]);
					if (fk.IsEqual(fj))
					{
						EFLi[j] = nfi[k];
						break;
					}
				}
			}
		}
		else return false;
	}
	return true;
}

//=============================================================================
FEFaceEdgeList::FEFaceEdgeList()
{

}

bool FEFaceEdgeList::Create(FEFaceList& faceList, FEEdgeList& edgeList)
{
	FENodeEdgeList NEL;
	NEL.Create(edgeList);

	int faces = faceList.Faces();
	m_FEL.resize(faces);
	for (int i = 0; i < faces; ++i)
	{
		vector<int>& edges = m_FEL[i];
		edges.clear();

		FEFaceList::FACE face = faceList.Face(i);

		// find the corresponding edges
		int n = face.ntype;
		for (int j = 0; j < n; ++j)
		{
			int a = face.node[j];
			int b = face.node[(j + 1) % n];

			int edge = -1;
			const std::vector<int>& a_edges = NEL.EdgeList(a);
			for (int k = 0; k < a_edges.size(); ++k)
			{
				const FEEdgeList::EDGE& ek = edgeList[a_edges[k]];
				if (((ek.node[0] == a) && (ek.node[1] == b)) ||
					((ek.node[1] == a) && (ek.node[0] == b)))
				{
					edge = a_edges[k];
					break;
				}
			}

			assert(edge >= 0);
			if (edge == -1) return false;
			edges.push_back(edge);
		}
	}

	return true;
}

int FEFaceEdgeList::Edges(int nface)
{
	return (int)m_FEL[nface].size();
}

const std::vector<int>& FEFaceEdgeList::EdgeList(int nface) const
{
	return m_FEL[nface];
}

//=============================================================================
FENodeEdgeList::FENodeEdgeList()
{

}

bool FENodeEdgeList::Create(FEEdgeList& edgeList)
{
	m_NEL.clear();
	FEMesh* mesh = edgeList.GetMesh();
	int N = mesh->Nodes();
	m_NEL.resize(N);

	int NE = edgeList.Edges();
	for (int i = 0; i < NE; ++i)
	{
		const FEEdgeList::EDGE& edge = edgeList[i];
		m_NEL[edge.node[0]].push_back(i);
		m_NEL[edge.node[1]].push_back(i);
	}

	return true;
}
