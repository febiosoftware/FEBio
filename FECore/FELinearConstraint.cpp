#include "stdafx.h"
#include "FELinearConstraint.h"
#include "FEMesh.h"
#include "FEModel.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
FELinearConstraint::FELinearConstraint(FEModel* pfem) : FEModelComponent(FEBC_ID, pfem) 
{
	m_off = 0.0;
}

//-----------------------------------------------------------------------------
FELinearConstraint::FELinearConstraint(const FELinearConstraint& LC) : FEModelComponent(FEBC_ID, LC.GetFEModel())
{
	master = LC.master;
	m_off = LC.m_off;
	int n = (int) LC.slave.size();
	vector<DOF>::const_iterator it = LC.slave.begin();
	for (int i=0; i<n; ++i, ++it) slave.push_back(*it);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::CopyFrom(const FELinearConstraint& LC)
{
	master = LC.master;
	m_off = LC.m_off;
	int n = (int)LC.slave.size();
	vector<DOF>::const_iterator it = LC.slave.begin();
	for (int i = 0; i<n; ++i, ++it) slave.push_back(*it);
}

//-----------------------------------------------------------------------------
void FELinearConstraint::SetMasterDOF(int dof, int node)
{
	master.dof = dof;
	master.node = node;
}

//-----------------------------------------------------------------------------
void FELinearConstraint::AddSlaveDof(int dof, int node, double v)
{
	DOF d;
	d.dof = dof;
	d.node = node;
	d.val = v;
	slave.push_back(d);
}

//-----------------------------------------------------------------------------
// Initialization.
// Make sure the master dof does not appear as a slave dof
bool FELinearConstraint::Init()
{
	int n = (int) slave.size();
	for (int i=0; i<n; ++i)
	{
		DOF& slaveNode = slave[i];
		if ((slaveNode.node == master.node) && (slaveNode.dof == master.dof)) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
// This is called during model activation (i.e. at the start of an analysis step)
// The master nodes are fixed in order to make sure that they are not assigned an equation number.
void FELinearConstraint::Activate()
{
	FEModelComponent::Activate();
	FEMesh& mesh = GetFEModel()->GetMesh();

	// we need the master node to be fixed so that no equation is allocated
	FENode& node = mesh.Node(master.node);
	node.m_BC[master.dof] = DOF_FIXED;
}

//-----------------------------------------------------------------------------
void FELinearConstraint::Deactivate()
{
	FEModelComponent::Deactivate();
	FEMesh& mesh = GetFEModel()->GetMesh();

	FENode& node = mesh.Node(master.node);
	node.m_BC[master.dof] = DOF_OPEN;
}

//-----------------------------------------------------------------------------
void FELinearConstraint::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		ar.write(&master, sizeof(DOF), 1);
		int n = (int) slave.size();
		ar << n;
		vector<DOF>::iterator it = slave.begin();
		for (int i=0; i<n; ++i, ++it) ar << it->val << it->node << it->dof;
	}
	else
	{
		slave.clear();
		ar.read(&master, sizeof(DOF), 1);
		int n;
		ar >> n;
		for (int i=0; i<n; ++i)
		{
			DOF dof;
			ar >> dof.val >> dof.node >> dof.dof;
			slave.push_back(dof);
		}
	}
}
