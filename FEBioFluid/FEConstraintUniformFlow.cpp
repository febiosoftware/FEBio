/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

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
#include "FEConstraintUniformFlow.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEMesh.h>

BEGIN_FECORE_CLASS(FEConstraintUniformFlow, FESurfaceConstraint)
    ADD_PARAMETER(m_lc.m_laugon, "laugon");
    ADD_PARAMETER(m_lc.m_tol, "tol");
    ADD_PARAMETER(m_lc.m_eps, "penalty");
    ADD_PARAMETER(m_lc.m_naugmin, "minaug");
    ADD_PARAMETER(m_lc.m_naugmax, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEConstraintUniformFlow::FEConstraintUniformFlow(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem), m_lc(pfem)
{
    m_vbar = 0;
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEConstraintUniformFlow::Activate()
{
    // don't forget to call base class
    FENLConstraint::Activate();
    m_lc.Activate();
}

//-----------------------------------------------------------------------------
bool FEConstraintUniformFlow::Init()
{
    // initialize surface
    m_surf.Init();
    
    // evaluate the nodal normals
    int N = m_surf.Nodes(), jp1, jm1;
    vec3d y[FEElement::MAX_NODES], n;
    m_nn.resize(N,vec3d(0,0,0));
    
    // loop over all elements
    for (int i=0; i<m_surf.Elements(); ++i)
    {
        FESurfaceElement& el = m_surf.Element(i);
        int ne = el.Nodes();
        
        // get the nodal coordinates
        for (int j=0; j<ne; ++j) y[j] = m_surf.Node(el.m_lnode[j]).m_rt;
        
        // calculate the normals
        for (int j=0; j<ne; ++j)
        {
            jp1 = (j+1)%ne;
            jm1 = (j+ne-1)%ne;
            n = (y[jp1] - y[j]) ^ (y[jm1] - y[j]);
            m_nn[el.m_lnode[j]] += n;
        }
    }
    
    // normalize all vectors
    for (int i=0; i<N; ++i) m_nn[i].unit();

    // get the dof indices
    int dof_wx = GetDOFIndex("wx");
    int dof_wy = GetDOFIndex("wy");
    int dof_wz = GetDOFIndex("wz");
    
    // create linear constraint
    // for a surface with prescribed velocity the constraint
    // on (vx, vy, vz) is
    //  nx*vx + ny*vy + nz*vz = vbar
    for (int i=0; i<N; ++i) {

        //  (1-nx^2)*vx - nx*ny*vy - nx*nz*vz = 0
        FEAugLagLinearConstraint* pLC = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        for (int j=0; j<3; ++j) {
            FENode& node = m_surf.Node(i);
            switch (j) {
            case 0: pLC->AddDOF(node.GetID(), dof_wx, m_nn[i].x); break;
            case 1: pLC->AddDOF(node.GetID(), dof_wy, m_nn[i].y); break;
            case 2: pLC->AddDOF(node.GetID(), dof_wz, m_nn[i].z); break;
            }
        }
        // add the linear constraint to the system
        m_lc.add(pLC);
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FEConstraintUniformFlow::Update()
{
    // evaluate the mean normal velocity on this surface

    // get the dof indices
    int dof_wx = GetDOFIndex("wx");
    int dof_wy = GetDOFIndex("wy");
    int dof_wz = GetDOFIndex("wz");

    /*
    m_vbar = 0;
    // loop over all elements
    for (int i=0; i<m_surf.Nodes(); ++i)
    {
        FENode& node = m_surf.Node(i);
        vec3d w = node.get_vec3d(dof_wx, dof_wy, dof_wz);
        double wn = w*m_nn[i];
        m_vbar += wn;
    }
    m_vbar /= m_surf.Nodes();    */

    m_vbar = 0;
    double area = 0;
    // loop over all elements
    for (int i=0; i<m_surf.Elements(); ++i)
    {
        FESurfaceElement& el = m_surf.Element(i);
        int ne = el.Nodes();
        vector<vec3d> w(ne,vec3d(0,0,0));
        
        // get the nodal fluid velocities
        for (int j=0; j<ne; ++j) w[j] = m_surf.Node(el.m_lnode[j]).get_vec3d(dof_wx, dof_wy, dof_wz);
        
        int nint = el.GaussPoints();
        double* gw = el.GaussWeights();
        vec3d t[2];
        for (int n=0; n<nint; ++n) {
            double* H = el.H(n);
            vec3d wn(0,0,0);
            for (int ie=0; ie<ne; ++ie) wn += w[ie]*H[ie];
            m_surf.CoBaseVectors(el, n, t);
            vec3d nu = t[0] ^ t[1];
            m_vbar += (wn*nu)*gw[n];
            area += nu.norm()*gw[n];
        }
    }
    
    m_vbar /= area;

    // assign mean velocity as rhs of constraint
    m_lc.m_rhs = m_vbar;
}

//-----------------------------------------------------------------------------
void FEConstraintUniformFlow::Serialize(DumpStream& ar) { m_lc.Serialize(ar); }
void FEConstraintUniformFlow::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FEConstraintUniformFlow::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FEConstraintUniformFlow::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FEConstraintUniformFlow::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
