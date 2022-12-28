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
#include "FEFixedBC.h"
#include "FENodeSet.h"
#include "FENode.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem) : FENodalBC(pfem)
{
}

FEFixedBC::FEFixedBC(FEModel* pfem, int dof, FENodeSet* ps) : FENodalBC(pfem)
{
	SetDOFList(dof);
	SetNodeSet(ps);
}

//-----------------------------------------------------------------------------
//! initialization
bool FEFixedBC::Init()
{
	if (GetNodeSet() == nullptr) return false;
	if (GetDofList().Size() == 0) return false;

	return FENodalBC::Init();
}

//-----------------------------------------------------------------------------
void FEFixedBC::Activate()
{
	FENodalBC::Activate();

	FENodeSet& nset = *GetNodeSet();
	int n = nset.Size();
	int dofs = m_dof.Size();
	for (int i = 0; i<n; ++i)
	{
		// make sure we only activate open dof's
		FENode& node = *nset.Node(i);
		for (int j=0; j<dofs; ++j)
		{
			int dofj = m_dof[j];
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
	FENodalBC::Deactivate();
	FENodeSet* pns = GetNodeSet();
	int N = pns->Size();
	size_t dofs = m_dof.Size();
	for (size_t i = 0; i<N; ++i)
	{
		// get the node
		FENode& node = *pns->Node(i);

		// set the dof to open
		for (size_t j = 0; j < dofs; ++j)
		{
			node.set_bc(m_dof[j], DOF_OPEN);
		}
	}
}

//-----------------------------------------------------------------------------
void FEFixedBC::CopyFrom(FEBoundaryCondition* bc)
{
	FEFixedBC* fbc = dynamic_cast<FEFixedBC*>(bc);
	m_dof = fbc->GetDofList();
}

//=============================================================================
BEGIN_FECORE_CLASS(FEFixedDOF, FEFixedBC)
	ADD_PARAMETER(m_dofs, "dofs", 0, "$(dof_list)");
END_FECORE_CLASS();

FEFixedDOF::FEFixedDOF(FEModel* fem) : FEFixedBC(fem)
{

}

void FEFixedDOF::SetDOFS(const std::vector<int>& dofs)
{
	m_dofs = dofs;
}

bool FEFixedDOF::Init()
{
	SetDOFList(m_dofs);
	return FEFixedBC::Init();
}
