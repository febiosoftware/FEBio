/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "FENode.h"
#include "DumpStream.h"

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
	m_ID.assign(n, -1);
	m_BC.assign(n, 0);
	m_val_t.assign(n, 0.0);
	m_val_p.assign(n, 0.0);
	m_Fr.assign(n, 0.0);
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
	m_d0 = n.m_d0;
    m_dt = n.m_dt;
    m_dp = n.m_dp;

	m_nID = n.m_nID;
	m_rid = n.m_rid;
	m_nstate = n.m_nstate;

	m_ID = n.m_ID;
	m_BC = n.m_BC;
	m_val_t = n.m_val_t;
	m_val_p = n.m_val_p;
	m_Fr = n.m_Fr;
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
	m_d0 = n.m_d0;
    m_dt = n.m_dt;
    m_dp = n.m_dp;

	m_nID = n.m_nID;
	m_rid = n.m_rid;
	m_nstate = n.m_nstate;

	m_ID = n.m_ID;
	m_BC = n.m_BC;
	m_val_t = n.m_val_t;
	m_val_p = n.m_val_p;
	m_Fr = n.m_Fr;

	return (*this);
}

//-----------------------------------------------------------------------------
// Serialize
void FENode::Serialize(DumpStream& ar)
{
	ar & m_nID;
	ar & m_rt & m_at;
	ar & m_rp & m_vp & m_ap;
	ar & m_Fr;
	ar & m_val_t & m_val_p;
    ar & m_dt & m_dp;
	if (ar.IsShallow() == false)
	{
		ar & m_nstate;
		ar & m_ID;
		ar & m_BC;
		ar & m_r0;
		ar & m_rid;
		ar & m_d0;
	}
}

//-----------------------------------------------------------------------------
//! Update nodal values, which copies the current values to the previous array
void FENode::UpdateValues()
{
	m_val_p = m_val_t;
}
