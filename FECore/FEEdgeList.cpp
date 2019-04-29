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



#include "stdafx.h"
#include "FEEdgeList.h"
#include "FEMesh.h"
#include "FENodeNodeList.h"
#include "FEDomain.h"
#include "FEElementList.h"
#include <set>
using namespace std;

FEEdgeList::FEEdgeList() : m_mesh(nullptr)
{

}

int FEEdgeList::Edges() const
{
	return (int) m_edgeList.size();
}

const FEEdgeList::EDGE& FEEdgeList::operator[] (int i)
{
	return m_edgeList[i];
}

const FEEdgeList::EDGE& FEEdgeList::Edge(int i) const
{
	return m_edgeList[i];
}

FEMesh* FEEdgeList::GetMesh()
{
	return m_mesh;
}

struct edge_less
{
	bool operator ()(const pair<int, int>& lhs, const pair<int, int>& rhs)
	{
		int na0 = lhs.first;
		int na1 = lhs.second;
		if (na0 > na1) { int tmp = na0; na0 = na1; na1 = tmp; }

		int nb0 = rhs.first;
		int nb1 = rhs.second;
		if (nb0 > nb1) { int tmp = nb0; nb0 = nb1; nb1 = tmp; }

		if (na0 != nb0) return na0 < nb0;
		else return na1 < nb1;
	}
};

bool FEEdgeList::Create(FEMesh* pmesh)
{
	if (pmesh == nullptr) return false;
	m_mesh = pmesh;
	FEMesh& mesh = *pmesh;

	// create the node-node list
	FEElementList elemList(mesh);

	set<pair<int, int>, edge_less> edgeSet;

	const int ETET[6][2] = { { 0, 1 },{ 1, 2 },{ 2, 0 },{ 0, 3 },{ 1, 3 },{ 2, 3 } };
	const int EHEX[12][2] = { { 0, 1 },{ 1, 2 },{ 2, 3 },{ 3, 0 },{ 4, 5 },{ 5, 6 },{ 6, 7 },{ 7, 4 },{ 0, 4 },{ 1, 5 },{ 2, 6 },{ 3, 7 } };

	for (FEElementList::iterator it = elemList.begin(); it != elemList.end(); ++it)
	{
		FEElement& el = *it;

		if (el.Shape() == ET_TET4)
		{
			for (int i = 0; i < 6; ++i)
			{
				pair<int, int> edge;
				edge.first = el.m_node[ETET[i][0]];
				edge.second = el.m_node[ETET[i][1]];

				edgeSet.insert(edge);
			}
		}
		else if (el.Shape() == ET_HEX8)
		{
			for (int i = 0; i < 12; ++i)
			{
				pair<int, int> edge;
				edge.first = el.m_node[EHEX[i][0]];
				edge.second = el.m_node[EHEX[i][1]];

				edgeSet.insert(edge);
			}
		}
		else return false;
	}

	// copy set into a vector
	size_t edges = edgeSet.size();
	m_edgeList.resize(edges);
	set< pair<int, int>, edge_less>::iterator edgeIter = edgeSet.begin();
	for (size_t i = 0; i < edges; ++i, ++edgeIter)
	{
		const pair<int, int>& edge = *edgeIter;
		EDGE& Edge = m_edgeList[i];
		Edge.ntype = 2;
		int n0 = edge.first;
		int n1 = edge.second;
		if (n0 > n1) { int tmp = n0; n0 = n1; n1 = tmp; }
		Edge.node[0] = n0;
		Edge.node[1] = n1;
	}

	return true;
}

int FEEdgeList::FindEdge(int a, int b)
{
	for (int i = 0; i < (int)m_edgeList.size(); ++i)
	{
		EDGE& edge = m_edgeList[i];
		if ((edge.node[0] == a) && (edge.node[1] == b)) return i;
		if ((edge.node[0] == b) && (edge.node[1] == a)) return i;
	}
	return -1;
}

//=============================================================================

FEElementEdgeList::FEElementEdgeList()
{

}

int FEElementEdgeList::Edges(int elem) const
{
	return (int)m_EEL[elem].size();
}

const std::vector<int>& FEElementEdgeList::EdgeList(int elem) const
{
	return m_EEL[elem];
}

// NOTE: This only works for TET4 and HEX8 elements!
bool FEElementEdgeList::Create(FEElementList& elemList, FEEdgeList& edgeList)
{
	FEMesh& mesh = *edgeList.GetMesh();

	const int ETET[6][2] = { { 0, 1 },{ 1, 2 },{ 2, 0 },{ 0, 3 },{ 1, 3 },{ 2, 3 } };
	const int EHEX[12][2] = { { 0, 1 },{ 1, 2 },{ 2, 3 },{ 3, 0 },{ 4, 5 },{ 5, 6 },{ 6, 7 },{ 7, 4 },{ 0, 4 },{ 1, 5 },{ 2, 6 },{ 3, 7 } };

	int NN = mesh.Nodes();
	vector<pair<int, int> > NI;
	NI.resize(NN);
	for (int i = 0; i<NN; ++i) NI[i].second = 0;
	for (int i = edgeList.Edges() - 1; i >= 0; --i)
	{
		const FEEdgeList::EDGE& et = edgeList.Edge(i);
		NI[et.node[0]].first = i;
		NI[et.node[0]].second++;
	}

	int NE = mesh.Elements();
	m_EEL.resize(NE);
	int i = 0;
	for (FEElementList::iterator it = elemList.begin(); it != elemList.end(); ++it, ++i)
	{
		const FEElement& el = *it;
		vector<int>& EELi = m_EEL[i];
		if (el.Shape() == FE_Element_Shape::ET_TET4)
		{
			EELi.resize(6);
			for (int j = 0; j<6; ++j)
			{
				int n0 = el.m_node[ETET[j][0]];
				int n1 = el.m_node[ETET[j][1]];

				if (n1 < n0) { int nt = n1; n1 = n0; n0 = nt; }

				int l0 = NI[n0].first;
				int ln = NI[n0].second;
				for (int l = 0; l<ln; ++l)
				{
					assert(edgeList[l0 + l].node[0] == n0);
					if (edgeList[l0 + l].node[1] == n1)
					{
						EELi[j] = l0 + l;
						break;
					}
				}
			}
		}
		else if (el.Shape() == FE_Element_Shape::ET_HEX8)
		{
			EELi.resize(12);
			for (int j = 0; j<12; ++j)
			{
				int n0 = el.m_node[EHEX[j][0]];
				int n1 = el.m_node[EHEX[j][1]];

				if (n1 < n0) { int nt = n1; n1 = n0; n0 = nt; }

				int l0 = NI[n0].first;
				int ln = NI[n0].second;
				for (int l = 0; l<ln; ++l)
				{
					assert(edgeList[l0 + l].node[0] == n0);
					if (edgeList[l0 + l].node[1] == n1)
					{
						EELi[j] = l0 + l;
						break;
					}
				}
			}
		}
	}
	return true;
}
