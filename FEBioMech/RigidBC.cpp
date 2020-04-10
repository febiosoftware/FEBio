/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "RigidBC.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include "FERigidBody.h"
#include <FECore/FEMaterial.h>
#include <FECore/FELoadCurve.h>
#include "FEMechModel.h"
#include "FERigidMaterial.h"
#include <FECore/DumpStream.h>

REGISTER_SUPER_CLASS(FERigidBC, FERIGIDBC_ID);

BEGIN_FECORE_CLASS(FERigidNodeSet, FEModelComponent)
	ADD_PARAMETER(m_rigidMat, "rb");
	ADD_PARAMETER(m_nshellBC, "clamp_shells");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidNodeSet::FERigidNodeSet(FEModel* pfem) : FERigidBC(pfem)
{
	m_rigidMat = -1;
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
	m_rigidMat = rs.m_rigidMat;
	m_node = rs.m_node;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::operator = (const FERigidNodeSet& rs)
{
	m_rigidMat = rs.m_rigidMat;
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
	// Make sure the rigid material exists
	FEModel& fem = *GetFEModel();

	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));
	if (pm == 0) return false;
	if (pm->GetRigidBodyID() < 0) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::Activate()
{
	FERigidBC::Activate();

	// get the rigid body's ID
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));
	int rid = pm->GetRigidBodyID(); assert(rid >= 0);

	// assign the rigid body ID
	FEModelComponent::Activate();
	FEMesh& mesh = fem.GetMesh();
	for (size_t i=0; i<m_node.size(); ++i)
	{
		FENode& node = mesh.Node(m_node[i]);
		if (m_nshellBC == CLAMPED_SHELL)
		{
			if (node.HasFlags(FENode::SHELL)) 
				node.SetFlags(FENode::RIGID_CLAMP);
		}
		node.m_rid = rid;
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
void FERigidNodeSet::SetRigidMaterialID(int rid)
{
	m_rigidMat = rid;
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
	ar & m_node & m_rigidMat & m_nshellBC;
}

//=============================================================================
BEGIN_FECORE_CLASS(FERigidBodyFixedBC, FERigidBC)
	ADD_PARAMETER(m_rigidMat, "rb");
	ADD_PARAMETER(m_dofs, "dofs", 0, "Rx\0Ry\0Rz\0Ru\0Rv\0Rw\0");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidBodyFixedBC::FERigidBodyFixedBC(FEModel* pfem) : FERigidBC(pfem)
{
	m_rigidMat = -1;
	m_rb = -1;
	m_binit = false;
}

//-----------------------------------------------------------------------------
bool FERigidBodyFixedBC::Init()
{
	// Make sure the rigid material ID is valid
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat -1));
	if (pm == 0) return false;

	// make sure we have a valid dof
	for (int i = 0; i < m_dofs.size(); ++i)
	{
		int dof_i = m_dofs[i];
		if ((dof_i < 0) || (dof_i >= 6)) return false;
	}

	m_binit = true;
	return FERigidBC::Init();
}

//-----------------------------------------------------------------------------
void FERigidBodyFixedBC::Activate()
{
	FERigidBC::Activate();

	// Get the Rigidbody ID
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1)); assert(pm);
	m_rb = pm->GetRigidBodyID(); assert(m_rb >= 0);

	if (m_binit == false) Init();
	if (m_binit)
	{
		FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
		FERigidBody& RB = *fem.GetRigidBody(m_rb);

		// we only fix the open dofs. If a user accidentally applied a fixed and prescribed
		// rigid degree of freedom, then we make sure the prescribed takes precedence.
		for (int i = 0; i < m_dofs.size(); ++i)
		{
			int dof_i = m_dofs[i];
			if (RB.m_BC[dof_i] == DOF_OPEN) RB.m_BC[dof_i] = DOF_FIXED;
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyFixedBC::Deactivate()
{
	if (m_binit)
	{
		FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
		FERigidBody& RB = *fem.GetRigidBody(m_rb);

		// Since fixed rigid dofs can be overwritten by prescribed dofs, 
		// we have to make sure that this dof is actually a fixed dof.
		for (int i = 0; i < m_dofs.size(); ++i)
		{
			int dof_i = m_dofs[i];
			if (RB.m_BC[dof_i] == DOF_FIXED) RB.m_BC[dof_i] = DOF_OPEN;
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyFixedBC::Serialize(DumpStream& ar)
{
	FEModelComponent::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dofs & m_rigidMat & m_binit & m_rb;
}

//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FERigidBodyDisplacement, FERigidBC)
	ADD_PARAMETER(m_rigidMat, "rb");
	ADD_PARAMETER(m_dof, "dof", 0, "Rx\0Ry\0Rz\0Ru\0Rv\0Rw\0");
	ADD_PARAMETER(m_val, "value");
	ADD_PARAMETER(m_brel, "relative")
END_FECORE_CLASS();

FERigidBodyDisplacement::FERigidBodyDisplacement(FEModel* pfem) : FERigidBC(pfem)
{
	m_rigidMat = -1;
	m_dof = -1;
	m_val = 0.0;
	m_ref= 0.0; 
	m_brel = false; 
	m_binit = false;
}

//-----------------------------------------------------------------------------
bool FERigidBodyDisplacement::Init()
{
	FEModel& fem = *GetFEModel();
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rigidMat - 1));
	if (pm == 0) return false;
	m_rigidMat = pm->GetRigidBodyID(); if (m_rigidMat < 0) return false;

	// make sure we have a valid dof
	if ((m_dof < 0)||(m_dof >=6)) return false;

	m_binit = true;
	return true;
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Activate()
{
	// don't forget to call the base class
	FEModelComponent::Activate();

	if (m_binit == false) Init();

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& RB = *fem.GetRigidBody(m_rigidMat);

	// set some stuff
	RB.m_pDC[m_dof] = this;

	// mark the dof as prescribed
	RB.m_BC[m_dof] = DOF_PRESCRIBED;

	// set the relative offset
	m_ref = 0.0;
	if (m_brel)
	{
		switch (m_dof)
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
		FERigidBody& RB = *fem.GetRigidBody(m_rigidMat);

		// turn off the prescribed displacement
		RB.m_pDC[m_dof] = 0;
		RB.m_BC[m_dof] = DOF_OPEN;
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Serialize(DumpStream& ar)
{
	FERigidBC::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dof & m_rigidMat & m_ref & m_binit & m_brel;
}

//-----------------------------------------------------------------------------
double FERigidBodyDisplacement::Value()
{
	return m_val + m_ref;
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::InitTimeStep()
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& RB = *fem.GetRigidBody(GetID());
	int I = GetBC();
	RB.m_dul[I] = Value() - RB.m_Ut[I];
}

//=============================================================================
FERigidIC::FERigidIC(FEModel* fem) : FERigidBC(fem)
{

}

//=============================================================================
BEGIN_FECORE_CLASS(FERigidBodyVelocity, FERigidBC)
	ADD_PARAMETER(m_rid, "rb");
	ADD_PARAMETER(m_vel, "value");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidBodyVelocity::FERigidBodyVelocity(FEModel* pfem) : FERigidIC(pfem) 
{
	m_rid = -1;
	m_vel = vec3d(0, 0, 0);
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
	FERigidBody& RB = *fem.GetRigidBody(m_rid);

	RB.m_vp = RB.m_vt = m_vel;
}

//=============================================================================
BEGIN_FECORE_CLASS(FERigidBodyAngularVelocity, FERigidBC)
	ADD_PARAMETER(m_rid, "rb");
	ADD_PARAMETER(m_w, "value");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidBodyAngularVelocity::FERigidBodyAngularVelocity(FEModel* pfem) : FERigidIC(pfem)
{
	m_rid = -1;
	m_w = vec3d(0, 0, 0);
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
	FERigidBody& RB = *fem.GetRigidBody(m_rid);
	RB.m_wp = RB.m_wt = m_w;
}
