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
#include "FEInitialCondition.h"
#include "FEModel.h"
#include "FEMesh.h"
#include "DumpStream.h"

REGISTER_SUPER_CLASS(FEInitialCondition, FEIC_ID);

FEInitialCondition::FEInitialCondition(FEModel* pfem) : FEModelComponent(pfem)
{
}

//-----------------------------------------------------------------------------
// set the nodeset for this component
void FEInitialCondition::AddNodes(const FENodeSet& nset)
{
	m_nodeSet = nset;
}

//-----------------------------------------------------------------------------
// set the list of degrees of freedom
void FEInitialCondition::SetDOFList(const std::vector<int>& dofList)
{
	m_dofs = dofList;
}

//-----------------------------------------------------------------------------
void FEInitialCondition::Activate()
{
	FEModelComponent::Activate();
	if (m_dofs.empty()) return;

	int dofs = (int)m_dofs.size();
	std::vector<double> val(dofs, 0.0);

	int N = (int)m_nodeSet.size();
	for (size_t i = 0; i<N; ++i)
	{
		FENode& node = *m_nodeSet.Node(i);

		// get the nodal values
		NodalValues(i, val);
		
		for (int j = 0; j < dofs; ++j)
		{
			node.set(m_dofs[j], val[j]);
		}
	}
}

//-----------------------------------------------------------------------------
// serialization
void FEInitialCondition::Serialize(DumpStream& ar)
{
	FEModelComponent::Serialize(ar);
	if (ar.IsShallow()) return;

	ar & m_dofs;
	ar & m_nodeSet;
}

//======================================================================================
BEGIN_FECORE_CLASS(FEInitialDOF, FEInitialCondition)
	ADD_PARAMETER(m_dof, "dof", 0, "@dof_list");
	ADD_PARAMETER(m_data, "value");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEInitialDOF::FEInitialDOF(FEModel* pfem) : FEInitialCondition(pfem), m_data(FE_DOUBLE)
{
	m_dof = -1;
}

//-----------------------------------------------------------------------------
void FEInitialDOF::SetDOF(int ndof) { m_dof = ndof; }

//-----------------------------------------------------------------------------
bool FEInitialDOF::Init()
{
	if (m_dof == -1) return false;
	std::vector<int> dofs;
	dofs.push_back(m_dof);
	SetDOFList(dofs);
	return FEInitialCondition::Init();
}

//-----------------------------------------------------------------------------
void FEInitialDOF::Serialize(DumpStream& ar)
{
	FEInitialCondition::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dof;
}

//-----------------------------------------------------------------------------
void FEInitialDOF::AddNodes(const FENodeSet& set)
{
	FEInitialCondition::AddNodes(set);
	int N = set.size();
	m_data.Create(N, 0.0);
}

//-----------------------------------------------------------------------------
void FEInitialDOF::Add(int node, double value)
{
	m_nodeSet.add(node);
	m_data.Add(value);
}

//-----------------------------------------------------------------------------
// return the values for node i
void FEInitialDOF::NodalValues(int inode, std::vector<double>& values)
{
	values[0] = m_data.getValue(inode);
}

//==============================================================================
FEInitialBCVec3D::FEInitialBCVec3D(FEModel* pfem) : FEInitialCondition(pfem) 
{ 
	m_dof[0] = m_dof[1] = m_dof[2] = -1; 
}

//-----------------------------------------------------------------------------
void FEInitialBCVec3D::SetDOF(int d0, int d1, int d2)
{ 
	m_dof[0] = d0; m_dof[1] = d1; m_dof[2] = d2; 
}

//-----------------------------------------------------------------------------
bool FEInitialBCVec3D::Init()
{
	if ((m_dof[0] == -1) || (m_dof[1] == -1) || (m_dof[2] == -1)) return false;
	std::vector<int> dofs;
	dofs.push_back(m_dof[0]);
	dofs.push_back(m_dof[1]);
	dofs.push_back(m_dof[2]);
	SetDOFList(dofs);
	return FEModelComponent::Init();
}

//-----------------------------------------------------------------------------
void FEInitialBCVec3D::Serialize(DumpStream& ar)
{
	FEInitialCondition::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dof;
	ar & m_data;
}

//-----------------------------------------------------------------------------
void FEInitialBCVec3D::Add(int nid, const vec3d& v)
{ 
	m_nodeSet.add(nid);
	m_data.push_back(v);
}

//-----------------------------------------------------------------------------
void FEInitialBCVec3D::AddNodes(const FENodeSet& nset)
{
	FEInitialCondition::AddNodes(nset);
	int nsize = nset.size();
	m_data.resize(nsize, vec3d(0, 0, 0));
}

//-----------------------------------------------------------------------------
void FEInitialBCVec3D::NodalValues(int inode, std::vector<double>& values)
{
	values[0] = m_data[inode].x;
	values[1] = m_data[inode].y;
	values[2] = m_data[inode].z;
}
