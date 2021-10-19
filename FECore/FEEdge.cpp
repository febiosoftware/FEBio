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
#include "FEEdge.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FEEdge::FEEdge(FEModel* fem) : FEMeshPartition(FE_DOMAIN_EDGE, fem)
{
}

//-----------------------------------------------------------------------------
FEEdge::~FEEdge()
{
}

//-----------------------------------------------------------------------------
FENodeList FEEdge::GetNodeList()
{
	FEMesh* pm = GetMesh();
	FENodeList set(pm);

	vector<int> tag(pm->Nodes(), 0);
	for (int i = 0; i<Elements(); ++i)
	{
		FELineElement& el = Element(i);
		int ne = el.Nodes();
		for (int j = 0; j<ne; ++j)
		{
			if (tag[el.m_node[j]] == 0)
			{
				set.Add(el.m_node[j]);
				tag[el.m_node[j]] = 1;
			}
		}
	}
	return set;
}

//-----------------------------------------------------------------------------
bool FEEdge::Init()
{
	// make sure that there is an edge defined
	if (Elements() == 0) return false;

	// get the mesh to which this edge belongs
	FEMesh& mesh = *GetMesh();

	// This array is used to keep tags on each node
	vector<int> tag; tag.assign(mesh.Nodes(), -1);

	// let's find all nodes the edge needs
	int nn = 0;
	int ne = Elements();
	for (int i=0; i<ne; ++i)
	{
		FELineElement& el = Element(i);
		el.SetLocalID(i);

		for (int j=0; j<el.Nodes(); ++j)
		{
			// get the global node number
			int m = el.m_node[j];
		
			// create a local node number
			if (tag[m] == -1) tag[m] = nn++;

			// set the local node number
			el.m_lnode[j] = tag[m];
		}
	}

	// allocate node index table
	m_Node.resize(nn);

	// fill the node index table
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		if (tag[i] >= 0)
		{
			m_Node[tag[i]] = i;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEEdge::Create(int nelems, int elemType)
{ 
	m_Elem.resize(nelems); 
	for (int i = 0; i < nelems; ++i)
	{
		FELineElement& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	if (elemType != -1)
		for (int i=0; i<nelems; ++i) m_Elem[i].SetType(elemType);
}
