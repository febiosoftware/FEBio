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
BEGIN_FECORE_CLASS(FEInitialNodeSet, FEInitialCondition)
	ADD_PROPERTY(m_nodeSet, "node_set", FEProperty::Reference);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEInitialNodeSet::FEInitialNodeSet(FEModel* fem) : FEInitialCondition(fem)
{
	m_nodeSet = nullptr;
}

//-----------------------------------------------------------------------------
// set the nodeset for this component
void FEInitialNodeSet::SetNodeSet(FENodeSet* nset)
{
	m_nodeSet = nset;
}

//-----------------------------------------------------------------------------
// get the node set
FENodeSet* FEInitialNodeSet::GetNodeSet()
{
	return m_nodeSet;
}

//-----------------------------------------------------------------------------
// set the list of degrees of freedom
void FEInitialNodeSet::SetDOFList(const std::vector<int>& dofList)
{
	m_dofs = dofList;
}

//-----------------------------------------------------------------------------
bool FEInitialNodeSet::Init()
{
	if (m_nodeSet == nullptr) return false;
	return FEInitialCondition::Init();
}

//-----------------------------------------------------------------------------
void FEInitialNodeSet::Activate()
{
	FEModelComponent::Activate();
	if (m_dofs.empty()) return;

	int dofs = (int)m_dofs.size();
	std::vector<double> val(dofs, 0.0);

	int N = (int)m_nodeSet->Size();
	for (int i = 0; i<N; ++i)
	{
		FENode& node = *m_nodeSet->Node(i);

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
void FEInitialNodeSet::Serialize(DumpStream& ar)
{
	FEModelComponent::Serialize(ar);
	if (ar.IsShallow()) return;

	ar & m_dofs;
	ar & m_nodeSet;
}

//======================================================================================
BEGIN_FECORE_CLASS(FEInitialDOF, FEInitialNodeSet)
	ADD_PARAMETER(m_dof, "dof", 0, "@dof_list");
	ADD_PARAMETER(m_data, "value");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEInitialDOF::FEInitialDOF(FEModel* pfem) : FEInitialNodeSet(pfem), m_data(FE_DOUBLE)
{
	m_dof = -1;
}

//-----------------------------------------------------------------------------
FEInitialDOF::FEInitialDOF(FEModel* fem, int ndof, FENodeSet* nset) : FEInitialNodeSet(fem), m_data(FE_DOUBLE)
{
	SetDOF(ndof);
	SetNodeSet(nset);
}

//-----------------------------------------------------------------------------
void FEInitialDOF::SetDOF(int ndof) { m_dof = ndof; }

//-----------------------------------------------------------------------------
bool FEInitialDOF::Init()
{
	if (FEInitialNodeSet::Init() == false) return false;
	if (m_dof == -1) return false;
	std::vector<int> dofs;
	dofs.push_back(m_dof);
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
void FEInitialDOF::NodalValues(int inode, std::vector<double>& values)
{
	values[0] = m_data;
}

//==============================================================================
FEInitialBCVec3D::FEInitialBCVec3D(FEModel* pfem) : FEInitialNodeSet(pfem)
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
	if (FEInitialNodeSet::Init() == false) return false;
	if ((m_dof[0] == -1) || (m_dof[1] == -1) || (m_dof[2] == -1)) return false;
	std::vector<int> dofs;
	dofs.push_back(m_dof[0]);
	dofs.push_back(m_dof[1]);
	dofs.push_back(m_dof[2]);
	SetDOFList(dofs);
	return true;
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
void FEInitialBCVec3D::NodalValues(int inode, std::vector<double>& values)
{
	values[0] = m_data.x;
	values[1] = m_data.y;
	values[2] = m_data.z;
}

//-----------------------------------------------------------------------------
void FEInitialBCVec3D::SetValue(const vec3d& v)
{
	m_data = v;
}
