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
#include "FEBoundaryCondition.h"
#include "FEFacetSet.h"
#include "FEModel.h"

REGISTER_SUPER_CLASS(FEBoundaryCondition, FEBC_ID);

//-----------------------------------------------------------------------------
FEBoundaryCondition::FEBoundaryCondition(FEModel* pfem) : FEModelComponent(pfem), m_nodeSet(&pfem->GetMesh())
{
}

//-----------------------------------------------------------------------------
FEBoundaryCondition::~FEBoundaryCondition()
{
}

//-----------------------------------------------------------------------------
// set the dof list
void FEBoundaryCondition::SetDOFList(const std::vector<int>& dofs)
{
	m_dofs = dofs;
}

//-----------------------------------------------------------------------------
// get the dof list
const std::vector<int> FEBoundaryCondition::GetDOFList()
{
	return m_dofs;
}
//-----------------------------------------------------------------------------
// get the node set
const FENodeSet& FEBoundaryCondition::GetNodeSet() const
{
	return m_nodeSet;
}

//-----------------------------------------------------------------------------
void FEBoundaryCondition::AddNode(int n)
{
	m_nodeSet.add(n);
}

//-----------------------------------------------------------------------------
// assign a node set to the prescribed BC
void FEBoundaryCondition::AddNodes(const FENodeSet& set)
{
	m_nodeSet.add(set);
}

//-----------------------------------------------------------------------------
// assign a surface to the BC
void FEBoundaryCondition::AddNodes(const FEFacetSet& surf)
{
	FENodeSet nset = surf.GetNodeSet();
	AddNodes(nset);
}

//-----------------------------------------------------------------------------
//! serialization
void FEBoundaryCondition::Serialize(DumpStream& ar)
{
	FEModelComponent::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_nodeSet;
	ar & m_dofs;
}

//-----------------------------------------------------------------------------
void FEBoundaryCondition::Deactivate()
{
	FEModelComponent::Deactivate();

	int N = m_nodeSet.size();
	size_t dofs = m_dofs.size();
	for (size_t i = 0; i<N; ++i)
	{
		// get the node
		FENode& node = *m_nodeSet.Node(i);

		// set the dof to open
		for (size_t j = 0; j < dofs; ++j)
		{
			node.set_bc(m_dofs[j], DOF_OPEN);
		}
	}
}

//-----------------------------------------------------------------------------
//! fill the prescribed values
void FEBoundaryCondition::PrepStep(std::vector<double>& u, bool brel)
{

}
