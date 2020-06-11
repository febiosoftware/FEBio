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
#include "FEInitialCondition.h"
#include "DumpStream.h"
#include "FENode.h"

REGISTER_SUPER_CLASS(FEInitialCondition, FEIC_ID);

FEInitialCondition::FEInitialCondition(FEModel* pfem) : FEModelComponent(pfem)
{
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FENodalIC, FEInitialCondition)
	ADD_PROPERTY(m_nodeSet, "node_set", FEProperty::Reference);
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
	FEModelComponent::Activate();
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
	FEModelComponent::Serialize(ar);
	if (ar.IsShallow()) return;

	ar & m_dofs;
	ar & m_nodeSet;
}

//======================================================================================
BEGIN_FECORE_CLASS(FEInitialDOF, FENodalIC)
	ADD_PARAMETER(m_dof, "dof", 0, "@dof_list");
	ADD_PARAMETER(m_data, "value");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEInitialDOF::FEInitialDOF(FEModel* pfem) : FENodalIC(pfem), m_data(FE_DOUBLE)
{
	m_dof = -1;
}

//-----------------------------------------------------------------------------
FEInitialDOF::FEInitialDOF(FEModel* fem, int ndof, FENodeSet* nset) : FENodalIC(fem), m_data(FE_DOUBLE)
{
	SetDOF(ndof);
	SetNodeSet(nset);
}

//-----------------------------------------------------------------------------
void FEInitialDOF::SetDOF(int ndof) { m_dof = ndof; }

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
void FEInitialDOF::GetNodalValues(int inode, std::vector<double>& values)
{
	values[0] = m_data;
}
