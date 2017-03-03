//
//  FEBackFlowStabilization.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/2/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FEBackFlowStabilization.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_PARAMETER_LIST(FEBackFlowStabilization, FESurfaceLoad)
ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
ADD_PARAMETER(m_rho, FE_PARAM_DOUBLE, "density");
END_PARAMETER_LIST()

//-----------------------------------------------------------------------------
//! constructor
FEBackFlowStabilization::FEBackFlowStabilization(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_beta = 1.0;
    m_rho = 1.0;
    
    // get the degrees of freedom
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEBackFlowStabilization::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

void FEBackFlowStabilization::ElementStiffness(FESurfaceElement& el, matrix& ke)
{
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    mat3dd I(1);
    
    // gauss weights
    double* w = el.GaussWeights();
    
    // nodal coordinates
    FEMesh& mesh = *m_psurf->GetMesh();
    vec3d rt[FEElement::MAX_NODES], vt[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j) {
        rt[j] = mesh.Node(el.m_node[j]).m_rt;
        vt[j] = mesh.Node(el.m_node[j]).get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
    }
    
    // repeat over integration points
    ke.zero();
    for (int k=0; k<nint; ++k)
    {
        double* N = el.H(k);
        double* Gr = el.Gr(k);
        double* Gs = el.Gs(k);
        
        vec3d v(0,0,0);
        vec3d dxr(0,0,0), dxs(0,0,0);
        for (int i=0; i<neln; ++i)
        {
            v += vt[i]*N[i];
            dxr += rt[i]*Gr[i];
            dxs += rt[i]*Gs[i];
        }
        
        vec3d n = dxr ^ dxs;
        double da = n.unit();
        
        double vn = v*n;
        if (vn < 0) {
            mat3d K = (I*vn + (v & n))*(m_beta*m_rho*da*w[k]);
            
            // calculate stiffness component
            for (int i=0; i<neln; ++i)
                for (int j=0; j<neln; ++j)
                {
                    mat3d Kab = K*(N[i]*N[j]);
                    ke[3*i  ][3*j  ] -= Kab(0,0); ke[3*i  ][3*j+1] -= Kab(0,1); ke[3*i  ][3*j+2] -= Kab(0,2);
                    ke[3*i+1][3*j  ] -= Kab(1,0); ke[3*i+1][3*j+1] -= Kab(1,1); ke[3*i+1][3*j+2] -= Kab(1,2);
                    ke[3*i+2][3*j  ] -= Kab(2,0); ke[3*i+2][3*j+1] -= Kab(2,1); ke[3*i+2][3*j+2] -= Kab(2,2);
                }
        }
    }
}

//-----------------------------------------------------------------------------
//! calculates the element force

void FEBackFlowStabilization::ElementForce(FESurfaceElement& el, vector<double>& fe)
{
    // nr integration points
    int nint = el.GaussPoints();
    
    // nr of element nodes
    int neln = el.Nodes();
    
    mat3dd I(1);
    
    // nodal coordinates
    FEMesh& mesh = *m_psurf->GetMesh();
    vec3d rt[FEElement::MAX_NODES], vt[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j) {
        rt[j] = mesh.Node(el.m_node[j]).m_rt;
        vt[j] = mesh.Node(el.m_node[j]).get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
    }
    
    // repeat over integration points
    zero(fe);
    double* w  = el.GaussWeights();
    for (int j=0; j<nint; ++j)
    {
        double* N  = el.H(j);
        double* Gr = el.Gr(j);
        double* Gs = el.Gs(j);
        
        // traction at integration points
        vec3d v(0,0,0);
        vec3d dxr(0,0,0), dxs(0,0,0);
        for (int i=0; i<neln; ++i)
        {
            v += vt[i]*N[i];
            dxr += rt[i]*Gr[i];
            dxs += rt[i]*Gs[i];
        }
        
        vec3d n = dxr ^ dxs;
        double da = n.unit();
        
        // force vector
        double vn = v*n;
        if (vn < 0) {
            vec3d f = v*vn*(m_beta*m_rho*da*w[j]);
            
            for (int i=0; i<neln; ++i)
            {
                fe[3*i  ] += N[i]*f.x;
                fe[3*i+1] += N[i]*f.y;
                fe[3*i+2] += N[i]*f.z;
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FEBackFlowStabilization::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEBackFlowStabilization::UnpackLM(FEElement& el, vector<int>& lm)
{
    FEMesh& mesh = *GetSurface().GetMesh();
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
void FEBackFlowStabilization::StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver)
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
        ElementStiffness(el, ke);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
void FEBackFlowStabilization::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
    vector<double> fe;
    vector<int> lm;
    
    FESurface& surf = GetSurface();
    int npr = surf.Elements();
    for (int i=0; i<npr; ++i)
    {
        FESurfaceElement& el = m_psurf->Element(i);
        
        // calculate nodal normal tractions
        int neln = el.Nodes();
        vector<double> tn(neln);
        
        int ndof = 3*neln;
        fe.resize(ndof);
        
        ElementForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // add element force vector to global force vector
        R.Assemble(el.m_node, lm, fe);
    }
}
