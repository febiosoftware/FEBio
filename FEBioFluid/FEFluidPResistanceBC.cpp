//
//  FEFluidPResistanceBC.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/16/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#include "FEFluidPResistanceBC.h"
#include "FEFluid.h"
#include "FEFluidP.h"
#include "FEFluidFSI.h"
#include "stdafx.h"
#include "FECore/FEModel.h"
#include "FECore/DOFS.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidPResistanceBC, FESurfaceLoad)
ADD_PARAMETER(m_R , "R");
ADD_PARAMETER(m_p0, "pressure_offset");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidPResistanceBC::FEFluidPResistanceBC(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_R = 0.0;
    m_alpha = 1.0;
    m_p0 = 0;
    
    m_dofWX = pfem->GetDOFIndex("wx");
    m_dofWY = pfem->GetDOFIndex("wy");
    m_dofWZ = pfem->GetDOFIndex("wz");
    m_dofEF  = pfem->GetDOFIndex("ef" );
    
    m_dofWXP = pfem->GetDOFIndex("wxp");
    m_dofWYP = pfem->GetDOFIndex("wyp");
    m_dofWZP = pfem->GetDOFIndex("wzp");
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidPResistanceBC::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidPResistanceBC::Init()
{
    FEModelComponent::Init();
    
    FESurface* ps = &GetSurface();
    ps->Init();
    
    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidPResistanceBC::Activate()
{
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // mark node as having prescribed DOF
        node.set_bc(m_dofEF, DOF_PRESCRIBED);
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEFluidPResistanceBC::Update()
{
    // evaluate the flow rate
    double Q = FlowRate();
    
    // calculate the resistance pressure
    double p = m_R*Q + m_p0;
    
    // prescribe this pressure at the nodes
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofEF] < -1)
        {
            FENode& node = ps->Node(i);
            // set node as having prescribed DOF
            node.set(m_dofEF, p);
        }
    }
}

//-----------------------------------------------------------------------------
//! evaluate the flow rate across this surface
double FEFluidPResistanceBC::FlowRate()
{
    double Q = 0;
    
    vec3d rt[FEElement::MAX_NODES];
    vec3d vt[FEElement::MAX_NODES];
    
    for (int iel=0; iel<m_psurf->Elements(); ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);
        
        // nr integration points
        int nint = el.GaussPoints();
        
        // nr of element nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (int i=0; i<neln; ++i) {
            FENode& node = m_psurf->GetMesh()->Node(el.m_node[i]);
            rt[i] = node.m_rt*m_alpha + node.m_rp*(1-m_alpha);
            vt[i] = node.get_vec3d(m_dofWX, m_dofWY, m_dofWZ)*m_alphaf + node.get_vec3d(m_dofWXP, m_dofWYP, m_dofWZP)*(1-m_alphaf);
        }
        
        double* Nr, *Ns;
        double* N;
        double* w  = el.GaussWeights();
        
        vec3d dxr, dxs, v;
        
        // repeat over integration points
        for (int n=0; n<nint; ++n)
        {
            N  = el.H(n);
            Nr = el.Gr(n);
            Ns = el.Gs(n);
            
            // calculate the velocity and tangent vectors at integration point
            dxr = dxs = v = vec3d(0,0,0);
            for (int i=0; i<neln; ++i)
            {
                v += vt[i]*N[i];
                dxr += rt[i]*Nr[i];
                dxs += rt[i]*Ns[i];
            }
            
            vec3d normal = dxr ^ dxs;
            double q = normal*v;
            Q += q*w[n];
        }
    }
    
    return Q;
}
