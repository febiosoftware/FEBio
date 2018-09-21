#include "stdafx.h"
#include "RigidBC.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include "FERigidBody.h"
#include "FERigidSystem.h"
#include <FECore/FEMaterial.h>
#include <FECore/LoadCurve.h>
#include "FEMechModel.h"

BEGIN_PARAMETER_LIST(FERigidNodeSet, FEBoundaryCondition)
	ADD_PARAMETER(m_nshellBC, FE_PARAM_INT, "clamp_shells");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FERigidNodeSet::FERigidNodeSet(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem)
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
FERigidNodeSet::FERigidNodeSet(const FERigidNodeSet& rs) : FEBoundaryCondition(FEBC_ID, rs.GetFEModel())
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

	FEMaterial* pm = fem.GetMaterial(GetRigidID());
	if (pm->IsRigid() == false) return false;
	if (pm->GetRigidBodyID() < 0) return false;

	// assign correct rigid body ID's to rigid nodes
	SetRigidID(pm->GetRigidBodyID());

	return true;
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::Activate()
{
	FEBoundaryCondition::Activate();
	FEMesh& mesh = GetFEModel()->GetMesh();
	for (size_t i=0; i<m_node.size(); ++i)
	{
		FENode& node = mesh.Node(m_node[i]);
		if (m_nshellBC == CLAMPED_SHELL)
		{
			if (node.HasFlags(FENode::SHELL)) 
				node.SetFlags(node.Flags() | FENode::RIGID_CLAMP);
		}
		node.m_rid = m_rid;
	}
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::SetNodeSet(FENodeSet& ns)
{
	int N = ns.size();
	m_node.resize(N);
	for (int i=0; i<N; ++i) m_node[i] = ns[i];
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::Deactivate()
{
	FEBoundaryCondition::Deactivate();
	FEMesh& mesh = GetFEModel()->GetMesh();
	for (size_t i=0; i<m_node.size(); ++i)
	{
		FENode& node = mesh.Node(m_node[i]);
		if (m_nshellBC == CLAMPED_SHELL)
		{
			if (node.HasFlags(FENode::SHELL))
				node.SetFlags(node.Flags() & ~FENode::RIGID_CLAMP);
		}
		node.m_rid = -1;
	}
}

//-----------------------------------------------------------------------------
void FERigidNodeSet::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_node << m_rid;
		ar << m_nshellBC;
	}
	else
	{
		ar >> m_node >> m_rid;		
		ar >> m_nshellBC;
	}
}

//-----------------------------------------------------------------------------
FERigidBodyFixedBC::FERigidBodyFixedBC(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem)
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
	FEMaterial* pm = fem.GetMaterial(id-1);
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
	if (ar.IsShallow()) return;

	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << bc << id << m_binit;
	}
	else
	{
		ar >> bc >> id >> m_binit;
	}
}

//-----------------------------------------------------------------------------
FERigidBodyDisplacement::FERigidBodyDisplacement(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem)
{
	id = -1;
	ref= 0.0; 
	brel = false; 
	m_binit = false;
}

//-----------------------------------------------------------------------------
bool FERigidBodyDisplacement::Init()
{
	FEModel& fem = *GetFEModel();
	FEMaterial* pm = fem.GetMaterial(id-1);
	id = pm->GetRigidBodyID(); if (id < 0) return false;

	// make sure we have a valid dof
	if ((bc < 0)||(bc>=6)) return false;

	m_binit = true;
	return true;
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Activate()
{
	// don't forget to call the base class
	FEBoundaryCondition::Activate();

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidSystem& rs = *fem.GetRigidSystem();
	FERigidBody& RB = *rs.Object(id);

	// set some stuff
	RB.m_pDC[bc] = this;

	// mark the dof as prescribed
	RB.m_BC[bc] = DOF_PRESCRIBED;

	// set the relative offset
	ref = 0.0;
	if (brel)
	{
		switch (bc)
		{
		case 0: ref = RB.m_rt.x - RB.m_r0.x; break;
		case 1: ref = RB.m_rt.y - RB.m_r0.y; break;
		case 2: ref = RB.m_rt.z - RB.m_r0.z; break;
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Deactivate()
{
	FEBoundaryCondition::Deactivate();

	// get the rigid body
	// Since Deactivate is called before Init (for multi-step analysis; in the FEBio input)
	// we have to make sure the data is initialized
	if (m_binit)
	{
		FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
		FERigidSystem& rs = *fem.GetRigidSystem();
		FERigidBody& RB = *rs.Object(id);

		// turn off the prescribed displacement
		RB.m_pDC[bc] = 0;
		RB.m_BC[bc] = DOF_OPEN;
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << bc << id << lc << sf << m_binit;
	}
	else
	{
		ar >> bc >> id >> lc >> sf >> m_binit;
	}
}

//-----------------------------------------------------------------------------
double FERigidBodyDisplacement::Value()
{
	FEModel& fem = *GetFEModel();
	if (lc < 0) return 0;
	else return sf*fem.GetLoadCurve(lc)->Value() + ref;
}

//-----------------------------------------------------------------------------
bool FERigidBodyVelocity::Init()
{
	FEModel& fem = *GetFEModel();
	FEMaterial* pm = fem.GetMaterial(m_rid-1);
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
	FEMaterial* pm = fem.GetMaterial(m_rid-1);
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
