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
#include "FEFluidRotationalVelocity.h"
#include "FECore/FEModel.h"
#include "FECore/FEElemElemList.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/log.h"
#include "FECore/LinearSolver.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidRotationalVelocity, FEPrescribedBC)
	ADD_PARAMETER(m_w, "angular_speed");
	ADD_PARAMETER(m_n, "axis"         );
	ADD_PARAMETER(m_p, "origin"       );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidRotationalVelocity::FEFluidRotationalVelocity(FEModel* pfem) : FEPrescribedBC(pfem)
{
    m_w = 0.0;
    m_n = vec3d(0,0,1);
    m_p = vec3d(0,0,0);
    
    m_dofWX = pfem->GetDOFIndex("wx");
    m_dofWY = pfem->GetDOFIndex("wy");
    m_dofWZ = pfem->GetDOFIndex("wz");
}

//-----------------------------------------------------------------------------
// copy data from another class
void FEFluidRotationalVelocity::CopyFrom(FEPrescribedBC* pbc)
{
	FEFluidRotationalVelocity* rv = dynamic_cast<FEFluidRotationalVelocity*>(pbc);
	m_w = rv->m_w;
	m_n = rv->m_n;
	m_p = rv->m_p;
	m_node = rv->m_node;
	CopyParameterListState(rv->GetParameterList());
}

//-----------------------------------------------------------------------------
//! add some nodes
void FEFluidRotationalVelocity::AddNodes(const FENodeSet& set)
{
	int N = set.size();
	m_node.resize(N);
	for (int i = 0; i<N; ++i) m_node[i] = set[i];
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidRotationalVelocity::Init()
{
    FEModelComponent::Init();
    
    m_n.unit();
    
    // evaluate nodal radial positions
	FEMesh& mesh = GetFEModel()->GetMesh();
    int N = m_node.size();
    m_r.resize(N,vec3d(0,0,0));
    for (int i=0; i<N; ++i) {
        vec3d x = mesh.Node(m_node[i]).m_r0 - m_p;
        m_r[i] = x - m_n*(x*m_n);
    }
    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidRotationalVelocity::Activate()
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	int N = m_node.size();
	for (int i=0; i<N; ++i)
    {
        FENode& node = mesh.Node(m_node[i]);
        // mark node as having prescribed DOF
        node.set_bc(m_dofWX, DOF_PRESCRIBED);
        node.set_bc(m_dofWY, DOF_PRESCRIBED);
        node.set_bc(m_dofWZ, DOF_PRESCRIBED);
    }
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidRotationalVelocity::Deactivate()
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	int N = m_node.size();
	for (int i = 0; i<N; ++i)
	{
		FENode& node = mesh.Node(m_node[i]);

		node.set_bc(m_dofWX, DOF_OPEN);
		node.set_bc(m_dofWY, DOF_OPEN);
		node.set_bc(m_dofWZ, DOF_OPEN);
	}
}

//-----------------------------------------------------------------------------
void FEFluidRotationalVelocity::PrepStep(std::vector<double>& ui, bool brel)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	int N = m_node.size(), I;
	for (int i = 0; i<N; ++i)
	{
		FENode& node = mesh.Node(m_node[i]);
		vec3d v = (m_n ^ m_r[i])*m_w;

		I = -node.m_ID[m_dofWX] - 2; if (I >= 0) ui[I] = (brel ? v.x - node.get(m_dofWX) : v.x);
		I = -node.m_ID[m_dofWY] - 2; if (I >= 0) ui[I] = (brel ? v.y - node.get(m_dofWY) : v.y);
		I = -node.m_ID[m_dofWZ] - 2; if (I >= 0) ui[I] = (brel ? v.z - node.get(m_dofWZ) : v.z);
	}
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the velocities
void FEFluidRotationalVelocity::Update()
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	int N = m_node.size();
	for (int i = 0; i<N; ++i)
	{
		FENode& node = mesh.Node(m_node[i]);
		// evaluate the velocity
        vec3d v = (m_n ^ m_r[i])*m_w;

		if (node.m_ID[m_dofWX] < -1) node.set(m_dofWX, v.x);
        if (node.m_ID[m_dofWY] < -1) node.set(m_dofWY, v.y);
        if (node.m_ID[m_dofWZ] < -1) node.set(m_dofWZ, v.z);
    }
}
