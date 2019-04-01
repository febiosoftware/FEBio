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
#include "FEFixedBC.h"
#include "FENodeSet.h"
#include "FEModel.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem) : FEBoundaryCondition(pfem)
{
	m_dof = -1;
}

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem, int node, int dof) : FEBoundaryCondition(pfem)
{
	m_node.push_back(node);
	m_dof = dof;
}

//-----------------------------------------------------------------------------
void FEFixedBC::AddNode(int node)
{
	m_node.push_back(node);
}

//-----------------------------------------------------------------------------
void FEFixedBC::AddNodes(const FENodeSet& ns)
{
	int N = ns.size();
	for (int i = 0; i<N; ++i) AddNode(ns[i]);
}

//-----------------------------------------------------------------------------
void FEFixedBC::SetDOF(int dof)
{
	m_dof = dof;
}

//-----------------------------------------------------------------------------
//! get the node list
std::vector<int> FEFixedBC::GetNodeList()
{
	return m_node;
}

//-----------------------------------------------------------------------------
void FEFixedBC::SetNodeList(const std::vector<int>& nodeList)
{
	m_node = nodeList;
}

//-----------------------------------------------------------------------------
void FEFixedBC::Serialize(DumpStream& ar)
{
	FEBoundaryCondition::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_node & m_dof;
}

//-----------------------------------------------------------------------------
void FEFixedBC::Activate()
{
	FEBoundaryCondition::Activate();
	if (m_dof >= 0)
	{
		FEMesh& mesh = GetFEModel()->GetMesh();
		int n = (int)m_node.size();
		for (int i = 0; i<n; ++i)
		{
			// make sure we only activate open dof's
			FENode& node = mesh.Node(m_node[i]);
			if (node.get_bc(m_dof) == DOF_OPEN)
			{
				node.set_bc(m_dof, DOF_FIXED);
				node.set(m_dof, 0.0);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFixedBC::Deactivate()
{
	FEBoundaryCondition::Deactivate();
	if (m_dof >= 0)
	{
		FEMesh& mesh = GetFEModel()->GetMesh();
		int n = (int)m_node.size();
		for (int i = 0; i<n; ++i)
		{
			FENode& node = mesh.Node(m_node[i]);
			node.set_bc(m_dof, DOF_OPEN);
		}
	}
}
