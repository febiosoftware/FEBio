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
#include "FEMechModel.h"
#include "FERigidBody.h"

BEGIN_FECORE_CLASS(FERigidNodeSet, FENodalBC)
	ADD_PARAMETER(m_rigidMat, "rb")->setEnums("$(rigid_materials)");
	ADD_PARAMETER(m_bshellBC, "clamp_shells")->SetFlags(FEParamFlag::FE_PARAM_HIDDEN);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidNodeSet::FERigidNodeSet(FEModel* pfem) : FENodalBC(pfem)
{
	m_rigidMat = -1;
	m_bshellBC = true;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::CopyFrom(FEBoundaryCondition* pbc)
{
	assert(false);
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

	return FENodalBC::Init();
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::Activate()
{
	FENodalBC::Activate();

	// get the rigid body's ID
	FEMechModel& fem = dynamic_cast<FEMechModel&>(*GetFEModel());
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));
	int rid = pm->GetRigidBodyID(); assert(rid >= 0);

	FERigidBody& rb = *fem.GetRigidBody(rid);

	// assign the rigid body ID
	FEStepComponent::Activate();
	FEMesh& mesh = fem.GetMesh();
	FENodeSet& nset = *GetNodeSet();
	for (size_t i=0; i<nset.Size(); ++i)
	{
		FENode& node = mesh.Node(nset[i]);
		if (m_bshellBC)
		{
			if (node.HasFlags(FENode::SHELL)) 
				node.SetFlags(FENode::RIGID_CLAMP);
		}
		node.m_rid = rid;

		// we need to take the current position of the node
		// and transform it into the reference configuration of the rigid body.
		vec3d dr = node.m_rt - rb.m_rt;
		quatd q = rb.GetRotation().Inverse();
		q.RotateVector(dr);
		node.m_ra = rb.m_r0 + dr;
	}
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::SetRigidMaterialID(int rid)
{
	m_rigidMat = rid;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::Deactivate()
{
	FENodalBC::Deactivate();
	FEMesh& mesh = GetFEModel()->GetMesh();
	FENodeSet& nset = *GetNodeSet();
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
	FENodalBC::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_rigidMat & m_bshellBC;
}
