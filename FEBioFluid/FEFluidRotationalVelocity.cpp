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
    
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
    m_dofE = pfem->GetDOFIndex("e");
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
        id = node.m_ID[m_dofVX]; if (id >= 0) node.m_ID[m_dofVX] = -id-2;
        id = node.m_ID[m_dofVY]; if (id >= 0) node.m_ID[m_dofVY] = -id-2;
        id = node.m_ID[m_dofVZ]; if (id >= 0) node.m_ID[m_dofVZ] = -id-2;
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
        if (node.m_ID[m_dofVX] < -1) node.set(m_dofVX, v.x);
        if (node.m_ID[m_dofVY] < -1) node.set(m_dofVY, v.y);
        if (node.m_ID[m_dofVZ] < -1) node.set(m_dofVZ, v.z);
    }
}
