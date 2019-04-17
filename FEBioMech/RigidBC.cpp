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
#include "RigidBC.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include "FERigidBody.h"
#include "FERigidSystem.h"
#include <FECore/FEMaterial.h>
#include <FECore/FELoadCurve.h>
#include "FEMechModel.h"
#include "FERigidMaterial.h"
#include <FECore/DumpStream.h>

REGISTER_SUPER_CLASS(FERigidBC, FERIGIDBC_ID);

BEGIN_FECORE_CLASS(FERigidNodeSet, FEModelComponent)
	ADD_PARAMETER(m_nshellBC, "clamp_shells");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidNodeSet::FERigidNodeSet(FEModel* pfem) : FERigidBC(pfem)
{
	m_rid = -1;
	m_nshellBC = CLAMPED_SHELL;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::SetShellBC(SHELL_BC bc)
{
	m_nshellBC = bc;
}

//-----------------------------------------------------------------------------
FERigidNodeSet::FERigidNodeSet(const FERigidNodeSet& rs) : FERigidBC(rs.GetFEModel())
{
	m_rid = rs.m_rid;
	m_node = rs.m_node;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::operator = (const FERigidNodeSet& rs)
{
	m_rid = rs.m_rid;
	m_node = rs.m_node;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::AddNode(int nid)
{
	m_node.push_back(nid);
}

//-----------------------------------------------------------------------------
bool FERigidNodeSet::Init()
{
	FEModel& fem = *GetFEModel();

	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(GetRigidID()));
	if (pm == 0) return false;
	if (pm->GetRigidBodyID() < 0) return false;

	// assign correct rigid body ID's to rigid nodes
	SetRigidID(pm->GetRigidBodyID());

	return true;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::Activate()
{
	FEModelComponent::Activate();
	FEMesh& mesh = GetFEModel()->GetMesh();
	for (size_t i=0; i<m_node.size(); ++i)
	{
		FENode& node = mesh.Node(m_node[i]);
		if (m_nshellBC == CLAMPED_SHELL)
		{
			if (node.HasFlags(FENode::SHELL)) 
				node.SetFlags(FENode::RIGID_CLAMP);
		}
		node.m_rid = m_rid;
	}
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::SetNodeSet(FENodeSet& ns)
{
	int N = ns.Size();
	m_node.resize(N);
	for (int i=0; i<N; ++i) m_node[i] = ns[i];
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::Deactivate()
{
	FEModelComponent::Deactivate();
	FEMesh& mesh = GetFEModel()->GetMesh();
	for (size_t i=0; i<m_node.size(); ++i)
	{
		FENode& node = mesh.Node(m_node[i]);
		if (m_nshellBC == CLAMPED_SHELL)
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
	ar & m_node & m_rid & m_nshellBC;
}

//-----------------------------------------------------------------------------
FERigidBodyFixedBC::FERigidBodyFixedBC(FEModel* pfem) : FERigidBC(pfem)
{
	id = -1;
	bc = -1;
	m_binit = false;
}

//-----------------------------------------------------------------------------
bool FERigidBodyFixedBC::Init()
{
	// At this point, the id variable points to the material.
	// We need to associate it with a rigid body.
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(id-1));
	if (pm == 0) return false;
	id = pm->GetRigidBodyID(); if (id < 0) return false;

	// make sure we have a valid dof
	if ((bc < 0)||(bc>=6)) return false;

	m_binit = true;
	return true;
}

//-----------------------------------------------------------------------------
void FERigidBodyFixedBC::Activate()
{
	if (m_binit)
	{
		FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
		FERigidSystem& rs = *fem.GetRigidSystem();
		FERigidBody& RB = *rs.Object(id);

		// we only fix the open dofs. If a user accidentally applied a fixed and prescribed
		// rigid degree of freedom, then we make sure the prescribed takes precedence.
		if (RB.m_BC[bc] == DOF_OPEN) RB.m_BC[bc] = DOF_FIXED;
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyFixedBC::Deactivate()
{
	if (m_binit)
	{
		FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
		FERigidSystem& rs = *fem.GetRigidSystem();
		FERigidBody& RB = *rs.Object(id);

		// Since fixed rigid dofs can be overwritten by prescribed dofs, 
		// we have to make sure that this dof is actually a fixed dof.
		if (RB.m_BC[bc] == DOF_FIXED) RB.m_BC[bc] = DOF_OPEN;
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyFixedBC::Serialize(DumpStream& ar)
{
	FEModelComponent::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & bc & id & m_binit;
}

//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FERigidBodyDisplacement, FEModelComponent)
	ADD_PARAMETER(m_val, "value");
END_FECORE_CLASS();

FERigidBodyDisplacement::FERigidBodyDisplacement(FEModel* pfem) : FERigidBC(pfem)
{
	m_id = -1;
	m_val = 0.0;
	m_ref= 0.0; 
	m_brel = false; 
	m_binit = false;
}

//-----------------------------------------------------------------------------
bool FERigidBodyDisplacement::Init()
{
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_id - 1));
	if (pm == 0) return false;
	m_id = pm->GetRigidBodyID(); if (m_id < 0) return false;

	// make sure we have a valid dof
	if ((m_bc < 0)||(m_bc>=6)) return false;

	m_binit = true;
	return true;
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Activate()
{
	// don't forget to call the base class
	FEModelComponent::Activate();

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidSystem& rs = *fem.GetRigidSystem();
	FERigidBody& RB = *rs.Object(m_id);

	// set some stuff
	RB.m_pDC[m_bc] = this;

	// mark the dof as prescribed
	RB.m_BC[m_bc] = DOF_PRESCRIBED;

	// set the relative offset
	m_ref = 0.0;
	if (m_brel)
	{
		switch (m_bc)
		{
		case 0: m_ref = RB.m_rt.x - RB.m_r0.x; break;
		case 1: m_ref = RB.m_rt.y - RB.m_r0.y; break;
		case 2: m_ref = RB.m_rt.z - RB.m_r0.z; break;
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Deactivate()
{
	FEModelComponent::Deactivate();

	// get the rigid body
	// Since Deactivate is called before Init (for multi-step analysis; in the FEBio input)
	// we have to make sure the data is initialized
	if (m_binit)
	{
		FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
		FERigidSystem& rs = *fem.GetRigidSystem();
		FERigidBody& RB = *rs.Object(m_id);

		// turn off the prescribed displacement
		RB.m_pDC[m_bc] = 0;
		RB.m_BC[m_bc] = DOF_OPEN;
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Serialize(DumpStream& ar)
{
	FEModelComponent::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_bc & m_id & m_val & m_ref & m_binit;
}

//-----------------------------------------------------------------------------
double FERigidBodyDisplacement::Value()
{
	return m_val + m_ref;
}

//-----------------------------------------------------------------------------
bool FERigidBodyVelocity::Init()
{
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rid - 1));
	if (pm == 0) return false;
	m_rid = pm->GetRigidBodyID(); if (m_rid < 0) return false;
	return true;
}

//-----------------------------------------------------------------------------
void FERigidBodyVelocity::Activate()
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidSystem& rs = *fem.GetRigidSystem();
	FERigidBody& RB = *rs.Object(m_rid);

	RB.m_vp = RB.m_vt = m_vel;
}

//-----------------------------------------------------------------------------
bool FERigidBodyAngularVelocity::Init()
{
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rid - 1));
	if (pm == 0) return false;
	m_rid = pm->GetRigidBodyID(); if (m_rid < 0) return false;
	return true;
}

//-----------------------------------------------------------------------------
void FERigidBodyAngularVelocity::Activate()
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidSystem& rs = *fem.GetRigidSystem();
	FERigidBody& RB = *rs.Object(m_rid);
	RB.m_wp = RB.m_wt = m_w;
}
