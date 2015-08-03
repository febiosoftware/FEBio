#include "stdafx.h"
#include "BC.h"
#include "FEModel.h"
#include "FESolver.h"
#include "FERigidBody.h"

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
bool FERigidBodyFixedBC::Init()
{
	FEModel& fem = *GetFEModel();
	FEMaterial* pm = fem.GetMaterial(id-1);
	id = pm->GetRigidBodyID(); if (id < 0) return false;
	return true;
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
bool FERigidBodyDisplacement::Init()
{
	FEModel& fem = *GetFEModel();
	FEMaterial* pm = fem.GetMaterial(id-1);
	id = pm->GetRigidBodyID(); if (id < 0) return false;
	return true;
}

//-----------------------------------------------------------------------------
/*
void FERigidBodyDisplacement::Activate()
{
	// don't forget to call the base class
	FEBoundaryCondition::Activate();

	// get the rigid body
	FEModel& fem = *GetFEModel();
	FERigidBody& RB = static_cast<FERigidBody&>(*fem.Object(id));

	// set some stuff
	RB.m_pDC[bc] = this;
	RB.m_BC[bc] = 1;

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
*/

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
bool FERigidBodyVelocity::Init()
{
	FEModel& fem = *GetFEModel();
	FEMaterial* pm = fem.GetMaterial(id-1);
	id = pm->GetRigidBodyID(); if (id < 0) return false;
	return true;
}

//-----------------------------------------------------------------------------
bool FERigidBodyAngularVelocity::Init()
{
	FEModel& fem = *GetFEModel();
	FEMaterial* pm = fem.GetMaterial(id-1);
	id = pm->GetRigidBodyID(); if (id < 0) return false;
	return true;
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

//-----------------------------------------------------------------------------
double FERigidBodyDisplacement::Value()
{
	FEModel& fem = *GetFEModel();
	if (lc < 0) return 0;
	else return sf*fem.GetLoadCurve(lc)->Value() + ref;
}
