#include "stdafx.h"
#include "FENode.h"

//=============================================================================
// FENode
//-----------------------------------------------------------------------------
FENode::FENode()
{
	// set the default state
	m_nstate = 0;

	// rigid body data
	m_rid = -1;

	// default ID
	m_nID = -1;
}

//-----------------------------------------------------------------------------
void FENode::SetDOFS(int n)
{
	// initialize dof stuff
	m_ID.assign(n, DOF_FIXED);
	m_BC.assign(n, DOF_OPEN);
	m_val.assign(n, 0.0);
}

//-----------------------------------------------------------------------------
FENode::FENode(const FENode& n)
{
	m_r0 = n.m_r0;
	m_rt = n.m_rt;
	m_at = n.m_at;
	m_rp = n.m_rp;
	m_vp = n.m_vp;
	m_ap = n.m_ap;
	m_Fr = n.m_Fr;
	m_d0 = n.m_d0;

	m_nID = n.m_nID;
	m_rid = n.m_rid;
	m_nstate = n.m_nstate;

	m_ID = n.m_ID;
	m_BC = n.m_BC;
	m_val = n.m_val;
}

//-----------------------------------------------------------------------------
FENode& FENode::operator = (const FENode& n)
{
	m_r0 = n.m_r0;
	m_rt = n.m_rt;
	m_at = n.m_at;
	m_rp = n.m_rp;
	m_vp = n.m_vp;
	m_ap = n.m_ap;
	m_Fr = n.m_Fr;
	m_d0 = n.m_d0;

	m_nID = n.m_nID;
	m_rid = n.m_rid;
	m_nstate = n.m_nstate;

	m_ID = n.m_ID;
	m_BC = n.m_BC;
	m_val = n.m_val;

	return (*this);
}
