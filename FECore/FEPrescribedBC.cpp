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
#include "FEPrescribedBC.h"
#include "FEFacetSet.h"
#include "FEModel.h"

BEGIN_FECORE_CLASS(FEPrescribedBC, FEBoundaryCondition)
	ADD_PARAMETER(m_brelative, "relative");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPrescribedBC::FEPrescribedBC(FEModel* pfem) : FEBoundaryCondition(pfem), m_nodeSet(&pfem->GetMesh())
{
	m_brelative = false;
}

//-----------------------------------------------------------------------------
// set the relative flag
void FEPrescribedBC::SetRelativeFlag(bool br)
{ 
	m_brelative = br; 
}

//-----------------------------------------------------------------------------
void FEPrescribedBC::SetDOF(int dof)
{
	std::vector<int> dofList;
	dofList.push_back(dof);
	SetDOFList(dofList);
}

//-----------------------------------------------------------------------------
// set the dof list
void FEPrescribedBC::SetDOFList(const std::vector<int>& dofs)
{
	m_dofs = dofs;
}

//-----------------------------------------------------------------------------
// get the dof list
const std::vector<int> FEPrescribedBC::GetDOFList()
{
	return m_dofs;
}

//-----------------------------------------------------------------------------
// get the node set
const FENodeSet& FEPrescribedBC::GetNodeSet()
{
	return m_nodeSet;
}

void FEPrescribedBC::AddNode(int n)
{
	m_nodeSet.add(n);
}

//-----------------------------------------------------------------------------
// assign a node set to the prescribed BC
void FEPrescribedBC::AddNodes(const FENodeSet& set)
{
	m_nodeSet.add(set);
}

//-----------------------------------------------------------------------------
// assign a surface to the BC
void FEPrescribedBC::AddNodes(const FEFacetSet& surf)
{
	FENodeSet nset = surf.GetNodeSet();
	AddNodes(nset);
}

//-----------------------------------------------------------------------------
void FEPrescribedBC::Activate()
{
	FEBoundaryCondition::Activate();

	int N = m_nodeSet.size();
	size_t dofs = m_dofs.size();
	if (m_brelative) m_rval.assign(N*dofs, 0.0);
	for (size_t i = 0; i<N; ++i)
	{
		// get the node
		FENode& node = *m_nodeSet.Node(i);

		// set the dofs to prescribed
		for (size_t j = 0; j < dofs; ++j)
		{
			node.set_bc(m_dofs[j], DOF_PRESCRIBED);

			if (m_brelative)
			{
				m_rval[i*dofs + j] = node.get(m_dofs[j]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
// return the value for node i, dof j
void FEPrescribedBC::NodalValues(int nodelid, std::vector<double>& val)
{
}

//-----------------------------------------------------------------------------
void FEPrescribedBC::Deactivate()
{
	FEBoundaryCondition::Deactivate();

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
// This function is called when the solver needs to know the 
// prescribed dof values. The brel flag indicates wheter the total 
// value is needed or the value with respect to the current nodal dof value
void FEPrescribedBC::PrepStep(std::vector<double>& ui, bool brel)
{
	int N = m_nodeSet.size();
	size_t dofs = m_dofs.size();
	vector<double> val(dofs, 0.0);
	for (size_t i = 0; i<N; ++i)
	{
		// get the node
		FENode& node = *m_nodeSet.Node(i);

		// get the values
		NodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			double uj = val[j];
			if (m_brelative)
			{
				uj += m_rval[i*dofs + j];
			}

			int I = -node.m_ID[m_dofs[j]] - 2; 
			if (I >= 0) ui[I] = (brel ? uj - node.get(m_dofs[j]) : uj);
		}
	}
}

//-----------------------------------------------------------------------------
// This is called during nodal update and should be used to enforce the 
// nodal degrees of freedoms
void FEPrescribedBC::Update()
{
	int N = m_nodeSet.size();
	size_t dofs = m_dofs.size();
	std::vector<double> val(dofs, 0.0);
	for (size_t i = 0; i<N; ++i)
	{
		// get the node
		FENode& node = *m_nodeSet.Node(i);

		// get the values
		NodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			double uj = val[j];
			if (m_brelative)
			{
				uj += m_rval[i*dofs + j];
			}

			node.set(m_dofs[j], uj);
		}
	}
}

//-----------------------------------------------------------------------------
// serialization
void FEPrescribedBC::Serialize(DumpStream& ar)
{
	FEBoundaryCondition::Serialize(ar);

	if (ar.IsShallow())
	{
		ar & m_rval;
	}
	else
	{
		ar & m_nodeSet;
		ar & m_dofs;
		ar & m_rval;
	}
}
