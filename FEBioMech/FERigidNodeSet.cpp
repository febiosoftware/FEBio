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
#include "FERigidNodeSet.h"
#include <FECore/FEModel.h>
#include "FERigidMaterial.h"

BEGIN_FECORE_CLASS(FERigidNodeSet, FEModelComponent)
	ADD_PARAMETER(m_rigidMat, "rb");
	ADD_PARAMETER(m_bshellBC, "clamp_shells");

	ADD_PROPERTY(m_nodeSet, "node_set", FEProperty::Reference);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidNodeSet::FERigidNodeSet(FEModel* pfem) : FEBoundaryCondition(pfem)
{
	m_rigidMat = -1;
	m_bshellBC = true;
	m_nodeSet = nullptr;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::CopyFrom(FEBoundaryCondition* pbc)
{
	assert(false);
}

//-----------------------------------------------------------------------------
FERigidNodeSet::FERigidNodeSet(const FERigidNodeSet& rs) : FEBoundaryCondition(rs.GetFEModel())
{
	m_rigidMat = rs.m_rigidMat;
	m_nodeSet = rs.m_nodeSet;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::operator = (const FERigidNodeSet& rs)
{
	m_rigidMat = rs.m_rigidMat;
	m_nodeSet = rs.m_nodeSet;
}

//-----------------------------------------------------------------------------
bool FERigidNodeSet::Init()
{
	// Make sure the rigid material exists
	FEModel& fem = *GetFEModel();
	if (m_rigidMat <= 0) return false;
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));
	if (pm == 0) return false;
	if (pm->GetRigidBodyID() < 0) return false;

	if (m_nodeSet == nullptr) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::Activate()
{
	FEBoundaryCondition::Activate();

	// get the rigid body's ID
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));
	int rid = pm->GetRigidBodyID(); assert(rid >= 0);

	// assign the rigid body ID
	FEModelComponent::Activate();
	FEMesh& mesh = fem.GetMesh();
	FENodeSet& nset = *m_nodeSet;
	for (size_t i=0; i<nset.Size(); ++i)
	{
		FENode& node = mesh.Node(nset[i]);
		if (m_bshellBC)
		{
			if (node.HasFlags(FENode::SHELL)) 
				node.SetFlags(FENode::RIGID_CLAMP);
		}
		node.m_rid = rid;
	}
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::SetNodeSet(FENodeSet* ns)
{
	m_nodeSet = ns;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::SetRigidMaterialID(int rid)
{
	m_rigidMat = rid;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::Deactivate()
{
	FEModelComponent::Deactivate();
	FEMesh& mesh = GetFEModel()->GetMesh();
	FENodeSet& nset = *m_nodeSet;
	for (size_t i=0; i<nset.Size(); ++i)
	{
		FENode& node = mesh.Node(nset[i]);
		if (m_bshellBC)
		{
			if (node.HasFlags(FENode::SHELL))
				node.UnsetFlags(FENode::RIGID_CLAMP);
		}
		node.m_rid = -1;
	}
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::Serialize(DumpStream& ar)
{
	FEModelComponent::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_nodeSet & m_rigidMat & m_bshellBC;
}
