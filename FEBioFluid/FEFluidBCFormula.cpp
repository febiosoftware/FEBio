//
//  FEFluidBCFormula.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/18/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#include "FEFluidBCFormula.h"
#include "FECore/FEModel.h"
#include "FECore/FEElemElemList.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/log.h"
#include "FECore/LinearSolver.h"

//=============================================================================
BEGIN_PARAMETER_LIST(FEFluidBCFormula, FESurfaceLoad)
    ADD_PARAMETER(m_s , FE_PARAM_DOUBLE, "scale"  );
    ADD_PARAMETER(m_sz, FE_PARAM_STRING, "formula");
    ADD_PARAMETER(m_bc, FE_PARAM_STRING, "bc"     );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEFluidBCFormula::FEFluidBCFormula(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_dofWX = pfem->GetDOFIndex("wx");
    m_dofWY = pfem->GetDOFIndex("wy");
    m_dofWZ = pfem->GetDOFIndex("wz");
    m_dofEF = pfem->GetDOFIndex("ef");
    
    m_dof = -1;
    m_s = 1;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidBCFormula::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidBCFormula::Init()
{
    FEModelComponent::Init();
    
    FESurface* ps = &GetSurface();
    ps->Init();
    
    if (strcmp(m_bc, "wx") == 0) m_dof = m_dofWX;
    else if (strcmp(m_bc, "wy") == 0) m_dof = m_dofWY;
    else if (strcmp(m_bc, "wz") == 0) m_dof = m_dofWZ;
    else if (strcmp(m_bc, "ef") == 0) m_dof = m_dofEF;
    
    return true;
}

//-----------------------------------------------------------------------------
//! Mark the nodes with prescribed boundary conditions
void FEFluidBCFormula::MarkBC()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();

    int id;
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        id = node.m_ID[m_dof]; if (id >= 0) node.m_ID[m_dof] = -id-2;
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the boundary conditions
void FEFluidBCFormula::SetBC()
{
    // prescribe this velocity at the nodes
    FESurface* ps = &GetSurface();
    
    int ierr;
    
    double t = GetFEModel()->GetTime().currentTime;
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // define a math parser object
        MathParser m;
        m.SetVariable("x", node.m_rt.x);
        m.SetVariable("y", node.m_rt.y);
        m.SetVariable("z", node.m_rt.z);
        m.SetVariable("X", node.m_r0.x);
        m.SetVariable("Y", node.m_r0.y);
        m.SetVariable("Z", node.m_r0.z);
        m.SetVariable("t", t);
        
        // evaluate the boundary condition
        double v = m.eval(m_sz, ierr)*m_s;
        if (node.m_ID[m_dof] < -1) node.set(m_dof, v);
    }
}
