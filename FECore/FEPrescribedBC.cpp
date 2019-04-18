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
FEPrescribedBC::FEPrescribedBC(FEModel* pfem) : FEBoundaryCondition(pfem), m_dofs(pfem)
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
// Set the node list
void FEPrescribedBC::SetNodeList(const FENodeList& nodeList)
{
	m_nodeList = nodeList;
}

//-----------------------------------------------------------------------------
void FEPrescribedBC::Activate()
{
	FEBoundaryCondition::Activate();

	int N = m_nodeList.Size();
	int dofs = m_dofs.Size();
	if (m_brelative) m_rval.assign(N*dofs, 0.0);
	for (int i = 0; i<N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

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
void FEPrescribedBC::Deactivate()
{
	FEBoundaryCondition::Deactivate();
	int N = m_nodeList.Size();
	int dofs = m_dofs.Size();
	for (int i = 0; i<N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

		// set the dof to open
		for (int j = 0; j < dofs; ++j)
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
	int N = m_nodeList.Size();
	int dofs = m_dofs.Size();
	vector<double> val(dofs, 0.0);
	for (int i = 0; i<N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

		// get the values
		GetNodalValues(i, val);
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
	int N = m_nodeList.Size();
	int dofs = m_dofs.Size();
	std::vector<double> val(dofs, 0.0);
	for (int i = 0; i<N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

		// get the values
		GetNodalValues(i, val);
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
// initialization
bool FEPrescribedBC::Init()
{
	// get the dof list from the derived class
	if (SetDofList(m_dofs) == false) return false;

	return FEBoundaryCondition::Init();
}

//-----------------------------------------------------------------------------
// serialization
void FEPrescribedBC::Serialize(DumpStream& ar)
{
	FEBoundaryCondition::Serialize(ar);
	ar & m_rval;
}

//=============================================================================

BEGIN_FECORE_CLASS(FEPrescribedNodeSet, FEPrescribedBC)
	ADD_PROPERTY(m_nodeSet, "node_set", FEProperty::Reference);
END_FECORE_CLASS();

FEPrescribedNodeSet::FEPrescribedNodeSet(FEModel* fem) : FEPrescribedBC(fem)
{
	m_nodeSet = nullptr;
}

void FEPrescribedNodeSet::SetNodeSet(FENodeSet* nodeSet)
{
	m_nodeSet = nodeSet;
}

const FENodeSet* FEPrescribedNodeSet::GetNodeSet()
{
	return m_nodeSet;
}

bool FEPrescribedNodeSet::Init()
{
	if (m_nodeSet == nullptr) return false;
	return FEPrescribedBC::Init();
}

void FEPrescribedNodeSet::Activate()
{
	if (m_nodeSet)
	{
		SetNodeList(m_nodeSet->GetNodeList());
	}
	FEPrescribedBC::Activate();
}

//=============================================================================

BEGIN_FECORE_CLASS(FEPrescribedSurface, FEPrescribedBC)
	ADD_PROPERTY(m_surface, "surface");
END_FECORE_CLASS();

FEPrescribedSurface::FEPrescribedSurface(FEModel* fem) : FEPrescribedBC(fem)
{
	m_surface = nullptr;
}

void FEPrescribedSurface::SetSurface(FEFacetSet* surface)
{
	m_surface = surface;
}

const FEFacetSet* FEPrescribedSurface::GetSurface()
{
	return m_surface;
}

bool FEPrescribedSurface::Init()
{
	if (m_surface == nullptr) return false;
	return FEPrescribedBC::Init();
}

void FEPrescribedSurface::Activate()
{
	if (m_surface)
	{
		SetNodeList(m_surface->GetNodeList());
	}
	FEPrescribedBC::Activate();
}
