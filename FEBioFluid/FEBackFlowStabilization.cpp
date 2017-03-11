//
//  FEBackFlowStabilization.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/2/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FEBackFlowStabilization.h"
#include "FECore/FEModel.h"
#include "FEFluid.h"

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_PARAMETER_LIST(FEBackFlowStabilization, FESurfaceLoad)
ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
ADD_PARAMETER(m_bout, FE_PARAM_BOOL, "outflow");
END_PARAMETER_LIST()

//-----------------------------------------------------------------------------
//! constructor
FEBackFlowStabilization::FEBackFlowStabilization(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_beta = 1.0;
    m_rho = 1.0;
    m_bout = true;
    m_dir = 1.0;
    m_k = 1.0;
    
    // get the degrees of freedom
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
    m_dofE  = pfem->GetDOFIndex("e" );
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEBackFlowStabilization::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEBackFlowStabilization::Init()
{
    FEModelComponent::Init();
    
    FESurface* ps = &GetSurface();
    ps->Init();
    // get fluid density from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = ps->Element(0);
    FEMesh* mesh = ps->GetMesh();
    FEElement* pe = mesh->FindElementFromID(el.m_elem[0]);
    if (pe == nullptr) return false;
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    FEFluid* fluid = dynamic_cast<FEFluid*> (pm);
    if (fluid == nullptr) return false;
    // get the density and bulk modulus
    m_rho = fluid->m_rhor;
    // get the bulk modulus
    FEMaterialPoint* fp = fluid->CreateMaterialPointData();
    m_k = fluid->BulkModulus(*fp);
    
    if (!m_bout) m_dir = -1;
    
    // evaluate surface normals
    vector<vec3d> sn(ps->Elements(),vec3d(0,0,0));
    m_nu.resize(ps->Nodes(),vec3d(0,0,0));
    vec3d r0[FEElement::MAX_NODES];
    for (int iel=0; iel<ps->Elements(); ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);
        
        // nr integration points
        int nint = el.GaussPoints();
        
        // nr of element nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (int i=0; i<neln; ++i)
            r0[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_r0;
        
        double* Nr, *Ns;
        double* w  = el.GaussWeights();
        
        vec3d dxr, dxs;
        
        // repeat over integration points
        for (int n=0; n<nint; ++n)
        {
            Nr = el.Gr(n);
            Ns = el.Gs(n);
            
            // calculate the tangent vectors at integration point
            dxr = dxs = vec3d(0,0,0);
            for (int i=0; i<neln; ++i)
            {
                dxr += r0[i]*Nr[i];
                dxs += r0[i]*Ns[i];
            }
            
            sn[iel] = dxr ^ dxs;
            sn[iel].unit();
        }
        
        // evaluate nodal normals by averaging surface normals
        for (int i=0; i<neln; ++i)
            m_nu[el.m_lnode[i]] += sn[iel];
    }

    for (int i=0; i<ps->Nodes(); ++i) m_nu[i].unit();
    
    return true;
}

//-----------------------------------------------------------------------------

void FEBackFlowStabilization::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
}

//-----------------------------------------------------------------------------
//! Mark the nodes with prescribed dilatation
void FEBackFlowStabilization::MarkDilatation()
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
//! Evaluate and prescribe the backflow pressure
void FEBackFlowStabilization::SetDilatation()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofE] < -1)
        {
            FENode& node = ps->Node(i);
            // extract the velocity
            vec3d v = node.get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
            double vn = v*m_nu[i];
            // evaluate dilatation
            double e = (vn*m_dir < 0) ? m_beta*m_rho/m_k*vn*vn*m_dir : 0;
            // set node as having prescribed DOF
            node.set(m_dofE, e);
        }
    }
}

