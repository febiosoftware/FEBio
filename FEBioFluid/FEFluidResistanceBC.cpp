//
//  FEFluidResistanceBC.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 9/28/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "FEFluidResistanceBC.h"
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
void FEFluidResistanceBC::UnpackLM(FEElement& el, vector<int>& lm)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    int N = el.Nodes();
    lm.resize(N*3);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = mesh.Node(n);
        vector<int>& id = node.m_ID;
        
        lm[3*i  ] = id[m_dofVX];
        lm[3*i+1] = id[m_dofVY];
        lm[3*i+2] = id[m_dofVZ];
    }
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
void FEFluidResistanceBC::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
    FEModel& fem = R.GetFEModel();
    
    vector<double> fe;
    vector<int> elm;
    
    vec3d r0[FEElement::MAX_NODES];
    
    // evaluate the flow rate
    double Q = FlowRate();
    
    // calculate the resistance pressure
    double p = m_R*Q;
    
    int i, n;
    int npr = (int)m_psurf->Elements();
    for (int iel=0; iel<npr; ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);
        
        int ndof = 3*el.Nodes();
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
            
            vec3d normal = dxr ^ dxs;
            vec3d f = normal*(-p*w[n]);
            
            for (i=0; i<neln; ++i)
            {
                fe[3*i  ] += N[i]*f.x;
                fe[3*i+1] += N[i]*f.y;
                fe[3*i+2] += N[i]*f.z;
            }
        }
        
        // get the element's LM vector and adjust it
        UnpackLM(el, elm);
        
        // add element force vector to global force vector
        R.Assemble(el.m_node, elm, fe);
    }
}

//-----------------------------------------------------------------------------
//! evaluate the stiffness matrix
void FEFluidResistanceBC::StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver)
{
    matrix ke;
    vector<int> lm;
    
    FESurface& surf = GetSurface();
    int npr = surf.Elements();
    for (int m=0; m<npr; ++m)
    {
        // get the surface element
        FESurfaceElement& el = m_psurf->Element(m);
        
        // calculate nodal normal tractions
        int neln = el.Nodes();
        vector<double> tn(neln);
        
        // get the element stiffness matrix
        int ndof = 3*neln;
        ke.resize(ndof, ndof);
        
        // calculate pressure stiffness
        FluidResistanceStiffness(el, ke);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! calculate stiffness for an element
void FEFluidResistanceBC::FluidResistanceStiffness(FESurfaceElement& el, matrix& ke)
{
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    // gauss weights
    double* w = el.GaussWeights();
    
    // nodal coordinates
    FEMesh& mesh = *m_psurf->GetMesh();
    vec3d rt[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;
    
    // repeat over integration points
    ke.zero();
    for (int n=0; n<nint; ++n)
    {
        double* N = el.H(n);
        double* Gr = el.Gr(n);
        double* Gs = el.Gs(n);
        
        vec3d dxr(0,0,0), dxs(0,0,0);
        for (int i=0; i<neln; ++i)
        {
            dxr += rt[i]*Gr[i];
            dxs += rt[i]*Gs[i];
        }
        
        vec3d nu = dxr ^ dxs;
        double da = nu.unit();
        
        // calculate stiffness component
        for (int i=0; i<neln; ++i)
            for (int j=0; j<neln; ++j)
            {
                mat3ds kab = dyad(nu)*(m_R*N[i]*N[j]*w[n]*da);
                
                ke.add(3*i, 3*j, kab);
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
            r0[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_r0;
            vt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
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
