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
