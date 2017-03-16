//
//  FEFluidNormalVelocity.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/2/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FEFluidNormalVelocity.h"
#include "FECore/FEModel.h"

//=============================================================================
BEGIN_PARAMETER_LIST(FEFluidNormalVelocity, FESurfaceLoad)
ADD_PARAMETER(m_velocity, FE_PARAM_DOUBLE    , "velocity");
ADD_PARAMETER(m_VC      , FE_PARAM_DATA_ARRAY, "value"   );
ADD_PARAMETER(m_bpv     , FE_PARAM_BOOL      , "prescribe_nodal_velocities");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEFluidNormalVelocity::FEFluidNormalVelocity(FEModel* pfem) : FESurfaceLoad(pfem), m_VC(FE_DOUBLE)
{
    m_velocity = 0.0;
    m_VC.set(1.0);
    m_bpv = true;
    
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
    m_dofE = pfem->GetDOFIndex("e");
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidNormalVelocity::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
    m_VC.Create(ps);
}

//-----------------------------------------------------------------------------
void FEFluidNormalVelocity::UnpackLM(FEElement& el, vector<int>& lm)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    int N = el.Nodes();
    lm.resize(N);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = mesh.Node(n);
        vector<int>& id = node.m_ID;
        
        lm[i] = id[m_dofE];
    }
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
void FEFluidNormalVelocity::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
    FEModel& fem = R.GetFEModel();
    
    vector<double> fe;
    vector<int> elm;
    
    vec3d r0[FEElement::MAX_NODES];
    
    int i, n;
    int npr = (int)m_VC.size();
    for (int iel=0; iel<npr; ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);
        
        int ndof = el.Nodes();
        fe.resize(ndof);
        
        // nr integration points
        int nint = el.GaussPoints();
        
        // nr of element nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (i=0; i<neln; ++i) r0[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_r0;
        
        double* Gr, *Gs;
        double* N;
        double* w  = el.GaussWeights();
        
        vec3d dxr, dxs;
        
        // get the velocity at the integration point
        double vn = m_VC.get<double>(iel)*m_velocity;
        
        // repeat over integration points
        zero(fe);
        for (n=0; n<nint; ++n)
        {
            N  = el.H(n);
            Gr = el.Gr(n);
            Gs = el.Gs(n);
            
            // calculate the tangent vectors
            dxr = dxs = vec3d(0,0,0);
            for (i=0; i<neln; ++i)
            {
                dxr.x += Gr[i]*r0[i].x;
                dxr.y += Gr[i]*r0[i].y;
                dxr.z += Gr[i]*r0[i].z;
                
                dxs.x += Gs[i]*r0[i].x;
                dxs.y += Gs[i]*r0[i].y;
                dxs.z += Gs[i]*r0[i].z;
            }
            
            double da = (dxr ^ dxs).norm();
            
            for (i=0; i<neln; ++i)
                fe[i] += N[i]*vn*w[n]*da;
        }
        
        // get the element's LM vector and adjust it
        UnpackLM(el, elm);
        
        // add element force vector to global force vector
        R.Assemble(el.m_node, elm, fe);
    }
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidNormalVelocity::Init()
{
    FEModelComponent::Init();
    
    FESurface* ps = &GetSurface();
    ps->Init();
    FEMesh* mesh = ps->GetMesh();

    // evaluate surface normals
    vector<vec3d> sn(ps->Elements(),vec3d(0,0,0));
    m_nu.resize(ps->Nodes(),vec3d(0,0,0));
    m_VN.resize(ps->Nodes(),0);
    vector<int> nf(ps->Nodes(),0);
    vec3d r0[FEElement::MAX_NODES];
    for (int iel=0; iel<ps->Elements(); ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);
        
        // nr integration points
        int nint = el.GaussPoints();
        
        // nr of element nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (int i=0; i<neln; ++i) {
            r0[i] = mesh->Node(el.m_node[i]).m_r0;
            m_VN[el.m_lnode[i]] += m_VC.get<double>(iel);
            ++nf[el.m_lnode[i]];
        }
        
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
    
    for (int i=0; i<ps->Nodes(); ++i) {
        m_nu[i].unit();
        m_VN[i] /= nf[i];
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! Mark the nodes with prescribed velocities
void FEFluidNormalVelocity::MarkVelocity()
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
void FEFluidNormalVelocity::SetVelocity()
{
    // prescribe this velocity at the nodes
    FESurface* ps = &GetSurface();
    
    int id;
    for (int i=0; i<ps->Nodes(); ++i)
    {
        // evaluate the velocity
        vec3d v = m_nu[i]*(m_velocity*m_VN[i]);
        FENode& node = ps->Node(i);
        if (node.m_ID[m_dofVX] < -1) node.set(m_dofVX, v.x);
        if (node.m_ID[m_dofVY] < -1) node.set(m_dofVY, v.y);
        if (node.m_ID[m_dofVZ] < -1) node.set(m_dofVZ, v.z);
    }
}
