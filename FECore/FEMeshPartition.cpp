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
#include "FEMeshPartition.h"
#include "FEMaterial.h"
#include "FEDataExport.h"
#include "FEMesh.h"
#include "DOFS.h"
#include "FEGlobalMatrix.h"
#include <string.h>
#include "FEModel.h"
#include "DumpStream.h"

REGISTER_SUPER_CLASS(FEMeshPartition, FEDOMAIN_ID);

//-----------------------------------------------------------------------------
FEMeshPartition::FEMeshPartition(int nclass, FEModel* fem) : FECoreBase(fem), m_pMesh(&fem->GetMesh()), m_nclass(nclass)
{
	m_bactive = true;
}

//-----------------------------------------------------------------------------
FEMeshPartition::~FEMeshPartition()
{
	// delete all data export classes
	if (m_Data.empty() == false)
	{
		size_t ND = m_Data.size();
		for (size_t i = 0; i<ND; ++i) delete m_Data[i];
		m_Data.clear();
	}
}

//-----------------------------------------------------------------------------
void FEMeshPartition::AddDataExport(FEDataExport* pd)
{
	if (pd) m_Data.push_back(pd);
}

//-----------------------------------------------------------------------------
FEElement* FEMeshPartition::FindElementFromID(int nid)
{
	for (int i = 0; i<Elements(); ++i)
	{
		FEElement& el = ElementRef(i);
		if (el.GetID() == nid) return &el;
	}

	return 0;
}

//-----------------------------------------------------------------------------
void FEMeshPartition::Serialize(DumpStream& ar)
{
	FECoreBase::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_Node;
	ar & m_dof;
	ar & m_nclass;
	ar & m_bactive;
//	ar & m_Data;
}

//-----------------------------------------------------------------------------
void FEMeshPartition::SetDOFList(vector<int>& dof)
{
	m_dof = dof;
}

//-----------------------------------------------------------------------------
// This is the default packing method. 
// It stores all the degrees of freedom for the first node in the order defined
// by the DOF array, then for the second node, and so on. 
void FEMeshPartition::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	int ndofs = (int)m_dof.size();
	lm.resize(N*ndofs);
	for (int i = 0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;
		for (int j = 0; j<ndofs; ++j) lm[i*ndofs + j] = id[m_dof[j]];
	}
}

//-----------------------------------------------------------------------------
void FEMeshPartition::BuildMatrixProfile(FEGlobalMatrix& M)
{
	vector<int> elm;
	const int NE = Elements();
	for (int j = 0; j<NE; ++j)
	{
		FEElement& el = ElementRef(j);
		UnpackLM(el, elm);
		M.build_add(elm);
	}
}

//-----------------------------------------------------------------------------
void FEMeshPartition::Activate()
{
	// get the number of degrees of freedom for this domain.
	const int ndofs = (int)m_dof.size();

	// activate all the degrees of freedom of this domain
	for (int i = 0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			for (int j = 0; j<ndofs; ++j)
				if (m_dof[j] >= 0) node.set_active(m_dof[j]);
		}
	}
}

//-----------------------------------------------------------------------------
//! return a specific node
FENode& FEMeshPartition::Node(int i)
{
	return m_pMesh->Node(m_Node[i]);
}

//-----------------------------------------------------------------------------
//! return a specific node
const FENode& FEMeshPartition::Node(int i) const
{
	return m_pMesh->Node(m_Node[i]);
}

//-----------------------------------------------------------------------------
void FEMeshPartition::CopyFrom(FEMeshPartition* pd)
{
	m_Node = pd->m_Node;
	m_dof = pd->m_dof;
	SetName(pd->GetName());
}

//-----------------------------------------------------------------------------
bool FEMeshPartition::Init()
{
	// base class first
	if (FECoreBase::Init() == false) return false;

	// make sure that there are elements in this domain
	if (Elements() == 0) return false;

	// get the mesh to which this domain belongs
	FEMesh& mesh = *GetMesh();

	// This array is used to keep tags on each node
	int NN = mesh.Nodes();
	vector<int> tag; tag.assign(NN, -1);

	// let's find all nodes the domain needs
	int nn = 0;
	int NE = Elements();
	for (int i = 0; i<NE; ++i)
	{
		FEElement& el = ElementRef(i);
		int ne = el.Nodes();
		for (int j = 0; j<ne; ++j)
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
	m_Node.assign(nn, -1);

	// fill the node index table
	for (int i = 0; i<NN; ++i)
	{
		if (tag[i] >= 0)
		{
			m_Node[tag[i]] = i;
		}
	}

#ifdef _DEBUG
	// make sure all nodes are assigned a local index
	for (int i = 0; i<nn; ++i)
	{
		assert(m_Node[i] >= 0);
	}
#endif

	return true;
}


//-----------------------------------------------------------------------------
void FEMeshPartition::ForEachMaterialPoint(std::function<void(FEMaterialPoint& mp)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEElement& el = ElementRef(i);
		int nint = el.GaussPoints();
		for (int n = 0; n < nint; ++n) f(*el.GetMaterialPoint(n));
	}
}

//-----------------------------------------------------------------------------
void FEMeshPartition::ForEachElement(std::function<void(FEElement& el)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i) f(ElementRef(i));
}
