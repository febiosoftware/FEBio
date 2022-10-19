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
#include "FEPrescribedBC.h"
#include "FESurface.h"
#include "FENode.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEPrescribedNodeSet, FENodalBC)
	ADD_PARAMETER(m_brelative, "relative");
END_FECORE_CLASS();

FEPrescribedNodeSet::FEPrescribedNodeSet(FEModel* fem) : FENodalBC(fem)
{
	m_brelative = false;
}

//-----------------------------------------------------------------------------
// set the relative flag
void FEPrescribedNodeSet::SetRelativeFlag(bool br)
{
	m_brelative = br;
}

void FEPrescribedNodeSet::Activate()
{
	FENodalBC::Activate();

	FENodeSet& nodeSet = *GetNodeSet();
	int N = nodeSet.Size();
	int dofs = m_dof.Size();
	if (m_brelative) m_rval.assign(N * dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *nodeSet.Node(i);

		// set the dofs to prescribed
		for (size_t j = 0; j < dofs; ++j)
		{
			node.set_bc(m_dof[j], DOF_PRESCRIBED);
			if (m_brelative)
			{
				m_rval[i * dofs + j] = node.get(m_dof[j]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEPrescribedNodeSet::Deactivate()
{
	FEBoundaryCondition::Deactivate();
	FENodeSet& nodeSet = *GetNodeSet();
	int N = nodeSet.Size();
	int dofs = m_dof.Size();
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *nodeSet.Node(i);

		// set the dof to open
		for (int j = 0; j < dofs; ++j)
		{
			node.set_bc(m_dof[j], DOF_OPEN);
		}
	}
}

//-----------------------------------------------------------------------------
// This function is called when the solver needs to know the 
// prescribed dof values. The brel flag indicates wheter the total 
// value is needed or the value with respect to the current nodal dof value
void FEPrescribedNodeSet::PrepStep(std::vector<double>& ui, bool brel)
{
	FENodeSet& nodeSet = *GetNodeSet();
	int N = nodeSet.Size();
	int dofs = m_dof.Size();
	vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *nodeSet.Node(i);

		// get the values
		GetNodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			double uj = val[j];
			if (m_brelative)
			{
				uj += m_rval[i * dofs + j];
			}

			int I = -node.m_ID[m_dof[j]] - 2;
			if (I >= 0) ui[I] = (brel ? uj - node.get(m_dof[j]) : uj);
		}
	}
}

//-----------------------------------------------------------------------------
// serialization
void FEPrescribedNodeSet::Serialize(DumpStream& ar)
{
	FENodalBC::Serialize(ar);
	ar & m_rval;
}

//-----------------------------------------------------------------------------
// This is called during nodal update and should be used to enforce the 
// nodal degrees of freedoms
void FEPrescribedNodeSet::Update()
{
	FENodeSet& nodeSet = *GetNodeSet();
	int N = nodeSet.Size();
	int dofs = m_dof.Size();
	std::vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *nodeSet.Node(i);

		// get the values
		GetNodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			double uj = val[j];
			if (m_brelative)
			{
				uj += m_rval[i * dofs + j];
			}

			node.set(m_dof[j], uj);
		}
	}
}

//-----------------------------------------------------------------------------
// This is called during contact update and should be used to enforce the
// nodal degrees of freedoms
void FEPrescribedNodeSet::Repair()
{
	FENodeSet& nodeSet = *GetNodeSet();
	int N = nodeSet.Size();
	int dofs = m_dof.Size();
	std::vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *nodeSet.Node(i);

		// get the values
		GetNodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			if (node.m_ID[m_dof[j]] >= 0) {
				node.m_ID[m_dof[j]] = -node.m_ID[m_dof[j]] - 2;
				double uj = val[j];
				if (m_brelative)
				{
					uj += m_rval[i * dofs + j];
				}

				node.set(m_dof[j], uj);
			}
		}
	}
}

//=============================================================================

BEGIN_FECORE_CLASS(FEPrescribedSurface, FESurfaceBC)
END_FECORE_CLASS();

FEPrescribedSurface::FEPrescribedSurface(FEModel* fem) : FESurfaceBC(fem)
{
	m_brelative = false;
}

void FEPrescribedSurface::Activate()
{
	FESurfaceBC::Activate();

	FESurface* surface = GetSurface();
	if (surface == nullptr) return;

	m_nodeList = surface->GetNodeList();
	
	int N = m_nodeList.Size();
	int dofs = m_dof.Size();
	if (m_brelative) m_rval.assign(N * dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

		// set the dofs to prescribed
		for (size_t j = 0; j < dofs; ++j)
		{
			node.set_bc(m_dof[j], DOF_PRESCRIBED);

			if (m_brelative)
			{
				m_rval[i * dofs + j] = node.get(m_dof[j]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEPrescribedSurface::Deactivate()
{
	FEBoundaryCondition::Deactivate();
	int N = m_nodeList.Size();
	int dofs = m_dof.Size();
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

		// set the dof to open
		for (int j = 0; j < dofs; ++j)
		{
			node.set_bc(m_dof[j], DOF_OPEN);
		}
	}
}

//-----------------------------------------------------------------------------
// This function is called when the solver needs to know the 
// prescribed dof values. The brel flag indicates wheter the total 
// value is needed or the value with respect to the current nodal dof value
void FEPrescribedSurface::PrepStep(std::vector<double>& ui, bool brel)
{
	int N = m_nodeList.Size();
	int dofs = m_dof.Size();
	vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
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
				uj += m_rval[i * dofs + j];
			}

			int I = -node.m_ID[m_dof[j]] - 2;
			if (I >= 0) ui[I] = (brel ? uj - node.get(m_dof[j]) : uj);
		}
	}
}

//-----------------------------------------------------------------------------
// set the relative flag
void FEPrescribedSurface::SetRelativeFlag(bool br)
{
	m_brelative = br;
}

//-----------------------------------------------------------------------------
// serialization
void FEPrescribedSurface::Serialize(DumpStream& ar)
{
	FEBoundaryCondition::Serialize(ar);
	ar & m_rval;
	if (ar.IsShallow() == false) ar & m_nodeList;
}

//-----------------------------------------------------------------------------
// This is called during nodal update and should be used to enforce the 
// nodal degrees of freedoms
void FEPrescribedSurface::Update()
{
	int N = m_nodeList.Size();
	int dofs = m_dof.Size();
	std::vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
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
				uj += m_rval[i * dofs + j];
			}

			node.set(m_dof[j], uj);
		}
	}
}

//-----------------------------------------------------------------------------
// This is called during contact update and should be used to enforce the
// nodal degrees of freedoms
void FEPrescribedSurface::Repair()
{
	int N = m_nodeList.Size();
	int dofs = m_dof.Size();
	std::vector<double> val(dofs, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the node
		FENode& node = *m_nodeList.Node(i);

		// get the values
		GetNodalValues(i, val);
		assert(val.size() == dofs);

		for (size_t j = 0; j < dofs; ++j)
		{
			if (node.m_ID[m_dof[j]] >= 0) {
				node.m_ID[m_dof[j]] = -node.m_ID[m_dof[j]] - 2;
				double uj = val[j];
				if (m_brelative)
				{
					uj += m_rval[i * dofs + j];
				}

				node.set(m_dof[j], uj);
			}
		}
	}
}
