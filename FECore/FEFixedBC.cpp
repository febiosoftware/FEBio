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

BEGIN_FECORE_CLASS(FEFixedBC, FEBoundaryCondition)
	ADD_PARAMETER(m_dofs, "dofs", 0, "@dof_list");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem) : FEBoundaryCondition(pfem)
{
}

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem, int node, int dof) : FEBoundaryCondition(pfem)
{
	AddNode(node);
	SetDOF(dof);
}

//-----------------------------------------------------------------------------
void FEFixedBC::Activate()
{
	FEBoundaryCondition::Activate();
	if (m_dofs.empty()) return;

	FENodeSet& nset = m_nodeSet;
	int n = nset.size();
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
void FEFixedBC::CopyFrom(FEBoundaryCondition* bc)
{
	FEFixedBC* fbc = dynamic_cast<FEFixedBC*>(bc);
	AddNodes(fbc->GetNodeSet());
	SetDOFList(fbc->GetDOFList());
}

//-----------------------------------------------------------------------------
void FEFixedBC::SetDOF(int ndof)
{
	std::vector<int> dofList;
	dofList.push_back(ndof);
	SetDOFList(dofList);
}
