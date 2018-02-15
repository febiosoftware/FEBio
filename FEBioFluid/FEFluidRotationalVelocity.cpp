//
//  FERotationalVelocity.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 6/2/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FEFluidRotationalVelocity.h"
#include "FECore/FEModel.h"
#include "FECore/FEElemElemList.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/log.h"
#include "FECore/LinearSolver.h"

//=============================================================================
BEGIN_PARAMETER_LIST(FEFluidRotationalVelocity, FESurfaceLoad)
ADD_PARAMETER(m_w       , FE_PARAM_DOUBLE    , "angular_speed");
ADD_PARAMETER(m_n       , FE_PARAM_VEC3D     , "axis"         );
ADD_PARAMETER(m_p       , FE_PARAM_VEC3D     , "origin"       );
END_PARAMETER_LIST();

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
//! Mark the nodes with prescribed velocities
void FEFluidRotationalVelocity::MarkVelocity()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();
    
    int id;
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        id = node.m_ID[m_dofWX]; if (id >= 0) node.m_ID[m_dofWX] = -id-2;
        id = node.m_ID[m_dofWY]; if (id >= 0) node.m_ID[m_dofWY] = -id-2;
        id = node.m_ID[m_dofWZ]; if (id >= 0) node.m_ID[m_dofWZ] = -id-2;
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the velocities
void FEFluidRotationalVelocity::SetVelocity()
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
