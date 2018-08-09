//
//  FEBackFlowResistance.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 8/8/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#include "FEBackFlowResistance.h"
#include "FEFluid.h"
#include "FEFluidFSI.h"
#include "stdafx.h"
#include "FECore/FEModel.h"
#include "FECore/DOFS.h"

//=============================================================================
BEGIN_PARAMETER_LIST(FEBackFlowResistance, FESurfaceLoad)
ADD_PARAMETER(m_R, FE_PARAM_DOUBLE      , "R");
ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE   , "beta");
ADD_PARAMETER(m_p0, FE_PARAM_DOUBLE     , "pressure_offset");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEBackFlowResistance::FEBackFlowResistance(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_beta = 0;
    m_R = 0;
    m_k = 1.0;
    m_rho = 0;
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
void FEBackFlowResistance::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEBackFlowResistance::Init()
{
    FEModelComponent::Init();
    
    FESurface* ps = &GetSurface();
    ps->Init();
    m_vn.resize(ps->Nodes());
    
    // get fluid bulk modulus from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = ps->Element(0);
    FEMesh* mesh = ps->GetMesh();
    FEElement* pe = mesh->FindElementFromID(el.m_elem[0]);
    if (pe == nullptr) return false;
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    FEFluid* fluid = dynamic_cast<FEFluid*> (pm);
    FEFluidFSI* fsi = dynamic_cast<FEFluidFSI*>(pm);
    // get the bulk modulus
    if (fluid) {
        FEMaterialPoint* fp = fluid->CreateMaterialPointData();
        m_k = fluid->BulkModulus(*fp);
        m_rho = fluid->m_rhor;
    }
    else if (fsi) {
        FEMaterialPoint* fp = fsi->CreateMaterialPointData();
        m_k = fsi->Fluid()->BulkModulus(*fp);
        m_rho = fsi->Fluid()->m_rhor;
    }
    else
        return false;
    
    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEBackFlowResistance::Activate()
{
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // mark node as having prescribed DOF
        node.m_BC[m_dofEF] = DOF_PRESCRIBED;
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEBackFlowResistance::Update()
{
    // evaluate the normal velocities
    NormalVelocities();
    
    // evaluate the flow rate
    double Q = FlowRate();
    
    // calculate the resistance pressure
    double pr = m_R*Q;
    
    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofEF] < -1)
        {
            // calculate the resistance pressure
            double p = (m_vn[i] < 0) ? -m_beta*m_rho*m_vn[i]*m_vn[i] : 0;

            // calculate the dilatation
            double e = -(pr+p+m_p0)/m_k;
            
            FENode& node = ps->Node(i);
            // set node as having prescribed DOF
            node.set(m_dofEF, e);
        }
    }
}

//-----------------------------------------------------------------------------
//! evaluate normal velocities
void FEBackFlowResistance::NormalVelocities()
{
    // evaluate node normals
    const int MN = FEElement::MAX_NODES;
    vec3d y[MN];
    int nnd = m_psurf->Nodes();
    vector<vec3d> nn(nnd,vec3d(0,0,0));

    // loop over all elements
    for (int i=0; i<m_psurf->Elements(); ++i)
    {
        FESurfaceElement& el = m_psurf->Element(i);
        int ne = el.Nodes();
        
        // get the nodal coordinates
        for (int j=0; j<ne; ++j) y[j] = m_psurf->GetMesh()->Node(el.m_lnode[j]).m_rt;
        
        // calculate the normals
        for (int j=0; j<ne; ++j)
        {
            int jp1 = (j+1)%ne;
            int jm1 = (j+ne-1)%ne;
            vec3d n = (y[jp1] - y[j]) ^ (y[jm1] - y[j]);
            nn[el.m_lnode[j]] += n;
        }
    }
    
    // normalize all vectors
    for (int i=0; i<nnd; ++i) nn[i].unit();
    
    // loop over all nodes and evaluate normal velocities
    for (int i=0; i<nnd; ++i) {
        FENode& node = m_psurf->GetMesh()->Node(i);
        vec3d vt = node.get_vec3d(m_dofWX, m_dofWY, m_dofWZ)*m_alphaf + node.get_vec3d(m_dofWXP, m_dofWYP, m_dofWZP)*(1-m_alphaf);
        m_vn[i] = vt*nn[i];
    }
    
}

//-----------------------------------------------------------------------------
//! evaluate the flow rate across this surface
double FEBackFlowResistance::FlowRate()
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
