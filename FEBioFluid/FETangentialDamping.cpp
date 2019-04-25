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
#include "FETangentialDamping.h"
#include <FECore/FEMesh.h>
#include "FEBioFluid.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FETangentialDamping, FESurfaceLoad)
	ADD_PARAMETER(m_eps, "penalty");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FETangentialDamping::FETangentialDamping(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
    m_eps = 0.0;
    
    // get the degrees of freedom
	m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETangentialDamping::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

void FETangentialDamping::ElementStiffness(FESurfaceElement& el, matrix& ke, const double alpha)
{
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    mat3dd I(1);
    
    // gauss weights
    double* w = el.GaussWeights();
    
    // nodal coordinates
    FEMesh& mesh = *m_psurf->GetMesh();
    vec3d rt[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;
    
    // repeat over integration points
    ke.zero();
    for (int k=0; k<nint; ++k)
    {
        double* N = el.H(k);
        double* Gr = el.Gr(k);
        double* Gs = el.Gs(k);
        
        vec3d dxr(0,0,0), dxs(0,0,0);
        for (int i=0; i<neln; ++i)
        {
            dxr += rt[i]*Gr[i];
            dxs += rt[i]*Gs[i];
        }
        
        vec3d n = dxr ^ dxs;
        double da = n.unit();
        
        mat3ds K = (I - dyad(n))*(-m_eps*da*w[k]);
        
        // calculate stiffness component
        for (int i=0; i<neln; ++i)
            for (int j=0; j<neln; ++j)
            {
                mat3ds Kab = K*(N[i]*N[j]);
                ke[3*i  ][3*j  ] += Kab.xx(); ke[3*i  ][3*j+1] += Kab.xy(); ke[3*i  ][3*j+2] += Kab.xz();
                ke[3*i+1][3*j  ] += Kab.xy(); ke[3*i+1][3*j+1] += Kab.yy(); ke[3*i+1][3*j+2] += Kab.yz();
                ke[3*i+2][3*j  ] += Kab.xz(); ke[3*i+2][3*j+1] += Kab.yz(); ke[3*i+2][3*j+2] += Kab.zz();
            }
    }
}

//-----------------------------------------------------------------------------
//! calculates the element force

void FETangentialDamping::ElementForce(FESurfaceElement& el, vector<double>& fe, const double alpha)
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
        
        // force vector
        vec3d f = (I - dyad(n))*v*(-m_eps*da*w[j]);
        
        for (int i=0; i<neln; ++i)
        {
            fe[3*i  ] += N[i]*f.x;
            fe[3*i+1] += N[i]*f.y;
            fe[3*i+2] += N[i]*f.z;
        }
    }
}

//-----------------------------------------------------------------------------
void FETangentialDamping::UnpackLM(FEElement& el, vector<int>& lm)
{
    FEMesh& mesh = *GetSurface().GetMesh();
    int N = el.Nodes();
    lm.resize(N*3);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = mesh.Node(n);
        vector<int>& id = node.m_ID;
        
        lm[3*i  ] = id[m_dofW[0]];
        lm[3*i+1] = id[m_dofW[1]];
        lm[3*i+2] = id[m_dofW[2]];
    }
}

//-----------------------------------------------------------------------------
void FETangentialDamping::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    FESurface& surf = GetSurface();
    int npr = surf.Elements();
#pragma omp parallel for
    for (int m=0; m<npr; ++m)
    {
		// get the surface element
		FESurfaceElement& el = m_psurf->Element(m);

		FEElementMatrix ke(el);
        
        // calculate nodal normal tractions
        int neln = el.Nodes();
        vector<double> tn(neln);
        
        // get the element stiffness matrix
        int ndof = 3*neln;
        ke.resize(ndof, ndof);
        
        // calculate pressure stiffness
        ElementStiffness(el, ke, tp.alpha);
        
        // get the element's LM vector
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
		LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
void FETangentialDamping::Residual(FEGlobalVector& R, const FETimeInfo& tp)
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
        
        int ndof = 3*neln;
        fe.resize(ndof);
        
        ElementForce(el, fe, tp.alpha);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // add element force vector to global force vector
        R.Assemble(el.m_node, lm, fe);
    }
}
