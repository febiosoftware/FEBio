#include "stdafx.h"
#include "BC.h"
#include "FEModel.h"
#include "FESolver.h"
#include "FERigidBody.h"
#include "log.h"

//-----------------------------------------------------------------------------
void FENodalForce::Serialize(DumpFile& ar)
{
	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << bc << lc << node << s;
	}
	else
	{
		ar >> bc >> lc >> node >> s;
	}
}

//-----------------------------------------------------------------------------
bool FENodalForce::Init()
{
	int NLC = GetFEModel()->LoadCurves();
	if ((lc < 0)||(lc >= NLC))
	{
		felog.printf("ERROR: Invalid loadcurve in nodal load %d\n", GetID());
		return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem)
{ 
	m_node = -1; 
	m_dof = -1; 
}

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem, int node, int dof) : FEBoundaryCondition(FEBC_ID, pfem)
{ 
	m_node = node; 
	m_dof = dof; 
}

//-----------------------------------------------------------------------------
void FEFixedBC::Serialize(DumpFile& ar)
{
	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_node << m_dof;
	}
	else
	{
		ar >> m_node >> m_dof;
	}
}

//-----------------------------------------------------------------------------
void FEFixedBC::Activate()
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	mesh.Node(m_node).m_ID[m_dof] = -1;
}

//-----------------------------------------------------------------------------
bool FEPrescribedBC::Init()
{
	// don't forget to call the base class
	if (FEBoundaryCondition::Init() == false) return false;

	// check the load curve ID
	FEModel& fem = *GetFEModel();
	int NLC = fem.LoadCurves();
	if ((lc < 0)||(lc >= NLC))
	{
		felog.printf("ERROR: Invalid loadcurve in prescribed BC %d\n", GetID());
		return false;
	}

	// make sure this is not a rigid node
	FEMesh& mesh = fem.GetMesh();
	if ((node < 0) || (node >= mesh.Nodes())) return false;
	if (mesh.Node(node).m_rid != -1) 
	{
		felog.printf("ERROR: Rigid nodes cannot be prescribed.\n");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEPrescribedBC::Serialize(DumpFile& ar)
{
	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << bc << lc << node << s << br << r;
	}
	else
	{
		ar >> bc >> lc >> node >> s >> br >> r;
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
		FEModel& fem = *GetFEModel();
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
		FEModel& fem = *GetFEModel();
		FERigidSystem& rs = *fem.GetRigidSystem();
		FERigidBody& RB = *rs.Object(id);

		// Since fixed rigid dofs can be overwritten by prescribed dofs, 
		// we have to make sure that this dof is actually a fixed dof.
		if (RB.m_BC[bc] == DOF_FIXED) RB.m_BC[bc] = DOF_OPEN;
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyFixedBC::Serialize(DumpFile& ar)
{
	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << bc << id;
	}
	else
	{
		ar >> bc >> id;
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
	FEModel& fem = *GetFEModel();
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
		FEModel& fem = *GetFEModel();
		FERigidSystem& rs = *fem.GetRigidSystem();
		FERigidBody& RB = *rs.Object(id);

		// turn off the prescribed displacement
		RB.m_pDC[bc] = 0;
		RB.m_BC[bc] = DOF_OPEN;
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyDisplacement::Serialize(DumpFile& ar)
{
	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << bc << id << lc << sf;
	}
	else
	{
		ar >> bc >> id >> lc >> sf;
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
	FEModel& fem = *GetFEModel();
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
	FEModel& fem = *GetFEModel();
	FERigidSystem& rs = *fem.GetRigidSystem();
	FERigidBody& RB = *rs.Object(m_rid);
	RB.m_wp = RB.m_wt = m_w;
}

//-----------------------------------------------------------------------------
void FERigidNode::Activate()
{
	FEBoundaryCondition::Activate();
	FEMesh& mesh = GetFEModel()->GetMesh();
	FENode& node = mesh.Node(nid);
	node.m_rid = rid;
}

//-----------------------------------------------------------------------------
void FERigidNode::Deactivate()
{
	FEBoundaryCondition::Deactivate();
	FEMesh& mesh = GetFEModel()->GetMesh();
	FENode& node = mesh.Node(nid);
	node.m_rid = -1;
}

//-----------------------------------------------------------------------------
void FERigidNode::Serialize(DumpFile& ar)
{
	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << nid << rid;		
	}
	else
	{
		ar >> nid >> rid;		
	}
}
