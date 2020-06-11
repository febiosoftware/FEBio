/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEFixedBC.h"
#include "FENodeSet.h"
#include "FEModel.h"
#include "DumpStream.h"

BEGIN_FECORE_CLASS(FEFixedBC, FEBoundaryCondition)
	ADD_PARAMETER(m_dofs, "dofs", 0, "@dof_list");
	ADD_PROPERTY(m_nodeSet, "node_set", FEProperty::Reference);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem) : FEBoundaryCondition(pfem)
{
	m_nodeSet = nullptr;
}

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem, int dof, FENodeSet* nset) : FEBoundaryCondition(pfem)
{
	SetDOF(dof);
	SetNodeSet(nset);
}

//-----------------------------------------------------------------------------
//! initialization
bool FEFixedBC::Init()
{
	if (m_nodeSet == nullptr) return false;
	if (m_dofs.size() == 0) return false;

	// set the DOF list
	m_dof.Clear();
	for (int i=0; i<(int)m_dofs.size(); ++i) m_dof.AddDof(m_dofs[i]);

	return FEBoundaryCondition::Init();
}

//-----------------------------------------------------------------------------
void FEFixedBC::Activate()
{
	FEBoundaryCondition::Activate();
	if (m_dofs.empty()) return;

	FENodeSet& nset = *m_nodeSet;
	int n = nset.Size();
	int dofs = m_dofs.size();
	for (int i = 0; i<n; ++i)
	{
		// make sure we only activate open dof's
		FENode& node = *nset.Node(i);
		for (int j=0; j<dofs; ++j)
		{
			int dofj = m_dofs[j];
			if (node.get_bc(dofj) == DOF_OPEN)
			{
				node.set_bc(dofj, DOF_FIXED);
				node.set(dofj, 0.0);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFixedBC::Deactivate()
{
	FEBoundaryCondition::Deactivate();
	int N = m_nodeSet->Size();
	size_t dofs = m_dofs.size();
	for (size_t i = 0; i<N; ++i)
	{
		// get the node
		FENode& node = *m_nodeSet->Node(i);

		// set the dof to open
		for (size_t j = 0; j < dofs; ++j)
		{
			node.set_bc(m_dofs[j], DOF_OPEN);
		}
	}
}

//-----------------------------------------------------------------------------
void FEFixedBC::CopyFrom(FEBoundaryCondition* bc)
{
	FEFixedBC* fbc = dynamic_cast<FEFixedBC*>(bc);
	SetDOFList(fbc->GetDOFList());
}

//-----------------------------------------------------------------------------
void FEFixedBC::SetDOF(int ndof)
{
	std::vector<int> dofList;
	dofList.push_back(ndof);
	SetDOFList(dofList);
}

//-----------------------------------------------------------------------------
void FEFixedBC::SetDOFList(const std::vector<int>& dofs)
{
	m_dofs = dofs;
}

//-----------------------------------------------------------------------------
// get the dof list
const std::vector<int> FEFixedBC::GetDOFList()
{
	return m_dofs;
}

//-----------------------------------------------------------------------------
void FEFixedBC::SetNodeSet(FENodeSet* nodeSet)
{
	m_nodeSet = nodeSet;
}

//-----------------------------------------------------------------------------
FENodeSet* FEFixedBC::GetNodeSet()
{
	return m_nodeSet;
}
