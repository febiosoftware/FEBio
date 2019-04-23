/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include "stdafx.h"
#include "FETangentialFlowStabilization.h"
#include "FEFluid.h"
#include "FEBioFluid.h"

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FETangentialFlowStabilization, FESurfaceLoad)
    ADD_PARAMETER(m_beta, "beta");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FETangentialFlowStabilization::FETangentialFlowStabilization(FEModel* pfem) : FESurfaceLoad(pfem), m_dofU(pfem), m_dofW(pfem)
{
    m_beta = 1.0;
    m_rho = 1.0;
    
    // get the degrees of freedom
	m_dofU.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::DISPLACEMENT));
	m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETangentialFlowStabilization::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FETangentialFlowStabilization::Init()
{
    FESurfaceLoad::Init();
    
    // get fluid density from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = m_psurf->Element(0);
    FEMesh* mesh = m_psurf->GetMesh();
    FEElement* pe = el.m_elem[0];
    if (pe == nullptr) return false;

	// get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    FEFluid* fluid = pm->ExtractProperty<FEFluid>();
	if (fluid == nullptr) return false;

	// get the density and bulk modulus
    m_rho = fluid->m_rhor;
    
    return true;
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

void FETangentialFlowStabilization::ElementStiffness(FESurfaceElement& el, matrix& ke, const double alpha)
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
        FENode& node = mesh.Node(el.m_node[j]);
        rt[j] = node.m_rt*alpha + node.m_rp*(1-alpha);
        vt[j] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*alpha + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1-alpha);
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
        
        vec3d vtau = (I - dyad(n))*v;
        double vmag = vtau.unit();
        mat3d K = (I - dyad(n) + dyad(vtau))*(-m_beta*m_rho*vmag*da*w[k]);
        // force vector (change sign for inflow vs outflow)
        vec3d ttt = vtau*(-m_beta*m_rho*vmag*w[k]);
        
        // calculate stiffness component
        for (int i=0; i<neln; ++i)
            for (int j=0; j<neln; ++j)
            {
                mat3d Kww = K*(N[i]*N[j]*alpha);
                vec3d g = (dxr*Gs[j] - dxs*Gr[j])*(N[i]*alpha);
                mat3d Kwu; Kwu.skew(g);
                Kwu = (ttt & n)*Kwu;
                ke[6*i+3][6*j  ] -= Kwu(0,0); ke[6*i+3][6*j+1] -= Kwu(0,1); ke[6*i+3][6*j+2] -= Kwu(0,2);
                ke[6*i+4][6*j  ] -= Kwu(1,0); ke[6*i+4][6*j+1] -= Kwu(1,1); ke[6*i+4][6*j+2] -= Kwu(1,2);
                ke[6*i+5][6*j  ] -= Kwu(2,0); ke[6*i+5][6*j+1] -= Kwu(2,1); ke[6*i+5][6*j+2] -= Kwu(2,2);
                ke[6*i+3][6*j+3] -= Kww(0,0); ke[6*i+3][6*j+4] -= Kww(0,1); ke[6*i+3][6*j+5] -= Kww(0,2);
                ke[6*i+4][6*j+3] -= Kww(1,0); ke[6*i+4][6*j+4] -= Kww(1,1); ke[6*i+4][6*j+5] -= Kww(1,2);
                ke[6*i+5][6*j+3] -= Kww(2,0); ke[6*i+5][6*j+4] -= Kww(2,1); ke[6*i+5][6*j+5] -= Kww(2,2);
            }
    }
}

//-----------------------------------------------------------------------------
//! calculates the element force

void FETangentialFlowStabilization::ElementForce(FESurfaceElement& el, vector<double>& fe, const double alpha)
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
        FENode& node = mesh.Node(el.m_node[j]);
        rt[j] = node.m_rt;
        vt[j] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*alpha + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1-alpha);
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
        
        // tangential traction = -beta*density*|tangential velocity|*(tangential velocity)
        vec3d vtau = (I - dyad(n))*v;
        double vmag = vtau.norm();
        // force vector (change sign for inflow vs outflow)
        vec3d f = vtau*(-m_beta*m_rho*vmag*da*w[j]);
        
        for (int i=0; i<neln; ++i)
        {
            fe[6*i+3] += N[i]*f.x;
            fe[6*i+4] += N[i]*f.y;
            fe[6*i+5] += N[i]*f.z;
        }
    }
}

//-----------------------------------------------------------------------------

void FETangentialFlowStabilization::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
	ar & m_rho;
}

//-----------------------------------------------------------------------------
void FETangentialFlowStabilization::UnpackLM(FEElement& el, vector<int>& lm)
{
    FEMesh& mesh = *GetSurface().GetMesh();
    int N = el.Nodes();
    lm.resize(N*6);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = mesh.Node(n);
        vector<int>& id = node.m_ID;
        
        lm[6*i  ] = id[m_dofU[0]];
        lm[6*i+1] = id[m_dofU[1]];
        lm[6*i+2] = id[m_dofU[2]];
        lm[6*i+3] = id[m_dofW[0]];
        lm[6*i+4] = id[m_dofW[1]];
        lm[6*i+5] = id[m_dofW[2]];
    }
}

//-----------------------------------------------------------------------------
void FETangentialFlowStabilization::StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp)
{
    FESurface& surf = GetSurface();
    int npr = surf.Elements();
#pragma omp parallel for
    for (int m=0; m<npr; ++m)
    {
        matrix ke;
        vector<int> lm;
        
        // get the surface element
        FESurfaceElement& el = m_psurf->Element(m);
        
        // calculate nodal normal tractions
        int neln = el.Nodes();
        vector<double> tn(neln);
        
        // get the element stiffness matrix
        int ndof = 6*neln;
        ke.resize(ndof, ndof);
        
        // calculate pressure stiffness
        ElementStiffness(el, ke, tp.alpha);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
void FETangentialFlowStabilization::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
    FESurface& surf = GetSurface();
    int npr = surf.Elements();
#pragma omp parallel for
    for (int i=0; i<npr; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        FESurfaceElement& el = m_psurf->Element(i);
        
        // calculate nodal normal tractions
        int neln = el.Nodes();
        vector<double> tn(neln);
        
        int ndof = 6*neln;
        fe.resize(ndof);
        
        ElementForce(el, fe, tp.alpha);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // add element force vector to global force vector
        R.Assemble(el.m_node, lm, fe);
    }
}
