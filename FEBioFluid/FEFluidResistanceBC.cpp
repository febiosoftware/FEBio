//
//  FEFluidResistanceBC.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 9/28/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "FEFluidResistanceBC.h"
#include "FEFluid.h"
#include "stdafx.h"
#include "FECore/FEModel.h"

//=============================================================================
BEGIN_PARAMETER_LIST(FEFluidResistanceBC, FESurfaceLoad)
ADD_PARAMETER(m_R, FE_PARAM_DOUBLE    , "R");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEFluidResistanceBC::FEFluidResistanceBC(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_R = 0.0;
    m_k = 1.0;
    m_alpha = 1.0;
    
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
    m_dofE  = pfem->GetDOFIndex("e" );
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidResistanceBC::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidResistanceBC::Init()
{
    FEModelComponent::Init();
    
    FESurface* ps = &GetSurface();
    ps->Init();
    // get fluid bulk modulus from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = ps->Element(0);
    FEMesh* mesh = ps->GetMesh();
    FEElement* pe = mesh->FindElementFromID(el.m_elem[0]);
    if (pe == nullptr) return false;
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    FEFluid* fluid = dynamic_cast<FEFluid*> (pm);
    if (fluid == nullptr) return false;
    // get the bulk modulus
    FEMaterialPoint* fp = fluid->CreateMaterialPointData();
    m_k = fluid->BulkModulus(*fp);
    
    return true;
}

//-----------------------------------------------------------------------------
//! Mark the nodes with prescribed dilatation
void FEFluidResistanceBC::MarkDilatation()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        int id = ps->Node(i).m_ID[m_dofE];
        if (id >= 0)
        {
            FENode& node = ps->Node(i);
            // mark node as having prescribed DOF
            node.m_ID[m_dofE] = -id-2;
        }
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEFluidResistanceBC::SetDilatation()
{
    // evaluate the flow rate
    double Q = FlowRate();
    
    // calculate the resistance pressure
    double p = m_R*Q;
    
    // calculate the dilatation
    double e = -p/m_k;
    
    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();

    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofE] < -1)
        {
            FENode& node = ps->Node(i);
            // set node as having prescribed DOF
            node.set(m_dofE, e);
        }
    }
}

//-----------------------------------------------------------------------------
//! evaluate the flow rate across this surface
double FEFluidResistanceBC::FlowRate()
{
    double Q = 0;
    
    vec3d r0[FEElement::MAX_NODES];
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
            r0[i] = node.m_r0;
            vt[i] = node.get_vec3d(m_dofVX, m_dofVY, m_dofVZ)*m_alpha + node.m_vp*(1-m_alpha);
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
                dxr += r0[i]*Nr[i];
                dxs += r0[i]*Ns[i];
            }
            
            vec3d normal = dxr ^ dxs;
            double q = normal*v;
            Q += q*w[n];
        }
    }

    return Q;
}
