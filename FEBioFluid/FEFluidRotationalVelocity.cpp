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
BEGIN_FECORE_CLASS(FEFluidRotationalVelocity, FESurfaceLoad)
	ADD_PARAMETER(m_w, "angular_speed");
	ADD_PARAMETER(m_n, "axis"         );
	ADD_PARAMETER(m_p, "origin"       );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidRotationalVelocity::FEFluidRotationalVelocity(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_w = 0.0;
    m_n = vec3d(0,0,1);
    m_p = vec3d(0,0,0);
    
    m_dofWX = pfem->GetDOFIndex("wx");
    m_dofWY = pfem->GetDOFIndex("wy");
    m_dofWZ = pfem->GetDOFIndex("wz");
    m_dofEF = pfem->GetDOFIndex("ef");
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidRotationalVelocity::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidRotationalVelocity::Init()
{
    FEModelComponent::Init();
    
    FESurface* ps = &GetSurface();
    ps->Init();
    
    m_n.unit();
    
    // evaluate nodal radial positions
    int N = ps->Nodes();
    m_r.resize(N,vec3d(0,0,0));
    for (int i=0; i<N; ++i) {
        vec3d x = ps->Node(i).m_r0 - m_p;
        m_r[i] = x - m_n*(x*m_n);
    }
    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidRotationalVelocity::Activate()
{
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // mark node as having prescribed DOF
        node.set_bc(m_dofWX, DOF_PRESCRIBED);
        node.set_bc(m_dofWY, DOF_PRESCRIBED);
        node.set_bc(m_dofWZ, DOF_PRESCRIBED);
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the velocities
void FEFluidRotationalVelocity::Update()
{
    // prescribe this velocity at the nodes
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        // evaluate the velocity
        vec3d v = (m_n ^ m_r[i])*m_w;
        FENode& node = ps->Node(i);
        if (node.m_ID[m_dofWX] < -1) node.set(m_dofWX, v.x);
        if (node.m_ID[m_dofWY] < -1) node.set(m_dofWY, v.y);
        if (node.m_ID[m_dofWZ] < -1) node.set(m_dofWZ, v.z);
    }
}
