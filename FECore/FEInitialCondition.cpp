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
#include "FEInitialCondition.h"
#include "DumpStream.h"
#include "FEMaterialPoint.h"
#include "FENode.h"
#include "FEModel.h"

FEInitialCondition::FEInitialCondition(FEModel* pfem) : FEStepComponent(pfem)
{
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FENodalIC, FEInitialCondition)
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENodalIC::FENodalIC(FEModel* fem) : FEInitialCondition(fem), m_dofs(fem)
{
	m_nodeSet = nullptr;
}

//-----------------------------------------------------------------------------
// set the nodeset for this component
void FENodalIC::SetNodeSet(FENodeSet* nset)
{
	m_nodeSet = nset;
}

//-----------------------------------------------------------------------------
// get the node set
FENodeSet* FENodalIC::GetNodeSet()
{
	return m_nodeSet;
}

//-----------------------------------------------------------------------------
// set the list of degrees of freedom
void FENodalIC::SetDOFList(const FEDofList& dofList)
{
	m_dofs = dofList;
}

//-----------------------------------------------------------------------------
bool FENodalIC::Init()
{
	if (m_nodeSet == nullptr) return false;
	return FEInitialCondition::Init();
}

//-----------------------------------------------------------------------------
void FENodalIC::Activate()
{
	FEStepComponent::Activate();
	if (m_dofs.IsEmpty()) return;

	int dofs = (int)m_dofs.Size();
	std::vector<double> val(dofs, 0.0);

	int N = (int)m_nodeSet->Size();
	for (int i = 0; i<N; ++i)
	{
		FENode& node = *m_nodeSet->Node(i);

		// get the nodal values
		GetNodalValues(i, val);
		
		for (int j = 0; j < dofs; ++j)
		{
			node.set(m_dofs[j], val[j]);
		}
	}
}

//-----------------------------------------------------------------------------
// serialization
void FENodalIC::Serialize(DumpStream& ar)
{
	FEStepComponent::Serialize(ar);
	if (ar.IsShallow()) return;

	ar & m_dofs;
	ar & m_nodeSet;
}

//======================================================================================
BEGIN_FECORE_CLASS(FEInitialDOF, FENodalIC)
	ADD_PARAMETER(m_dof, "dof", 0, "$(dof_list)");
	ADD_PARAMETER(m_data, "value");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEInitialDOF::FEInitialDOF(FEModel* pfem) : FENodalIC(pfem)
{
	m_dof = -1;
	m_data = 0.0;
}

//-----------------------------------------------------------------------------
FEInitialDOF::FEInitialDOF(FEModel* fem, int ndof, FENodeSet* nset) : FENodalIC(fem)
{
	SetDOF(ndof);
	SetNodeSet(nset);
	m_data = 0.0;
}

//-----------------------------------------------------------------------------
void FEInitialDOF::SetDOF(int ndof) { m_dof = ndof; }

//-----------------------------------------------------------------------------
bool FEInitialDOF::SetDOF(const char* szdof)
{
	FEModel* fem = GetFEModel();
	int ndof = fem->GetDOFIndex(szdof);
	assert(ndof >= 0);
	if (ndof < 0) return false;
	SetDOF(ndof);
	return true;
}

//-----------------------------------------------------------------------------
bool FEInitialDOF::Init()
{
	if (FENodalIC::Init() == false) return false;
	if (m_dof == -1) return false;
	FEDofList dofs(GetFEModel());
	if (dofs.AddDof(m_dof) == false) return false;
	SetDOFList(dofs);
	return true;
}

//-----------------------------------------------------------------------------
void FEInitialDOF::Serialize(DumpStream& ar)
{
	FEInitialCondition::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dof;
}

//-----------------------------------------------------------------------------
void FEInitialDOF::SetValue(double v)
{
	m_data = v;
}

//-----------------------------------------------------------------------------
// return the values for node i
void FEInitialDOF::GetNodalValues(int inode, std::vector<double>& val)
{
	assert(val.size() == 1);
	const FENodeSet& nset = *GetNodeSet();
	int nid = nset[inode];
	const FENode& node = *nset.Node(inode);

	FEMaterialPoint mp;
	mp.m_r0 = node.m_r0;
	mp.m_index = inode;

	val[0] = m_data(mp);
}
