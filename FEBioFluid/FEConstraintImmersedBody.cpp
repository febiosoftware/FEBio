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
#include "FEConstraintImmersedBody.h"
#include "FEBioFluid.h"
#include <FECore/FEElementList.h>
#include <FECore/FEEdgeList.h>

BEGIN_FECORE_CLASS(FEConstraintImmersedBody, FESurfaceConstraint)
//    ADD_PARAMETER(m_lc.m_laugon, "laugon");
//    ADD_PARAMETER(m_lc.m_tol, "tol");
//    ADD_PARAMETER(m_lc.m_eps, "penalty");
//    ADD_PARAMETER(m_lc.m_rhs, "rhs");
//    ADD_PARAMETER(m_lc.m_naugmin, "minaug");
//    ADD_PARAMETER(m_lc.m_naugmax, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEConstraintImmersedBody::FEConstraintImmersedBody(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem), m_lc(pfem), m_dofW(pfem), m_dofEF(pfem)
{
    // TODO: Can this be done in Init, since  there is no error checking
    if (pfem)
    {
        m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
        m_dofEF.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION));
    }
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEConstraintImmersedBody::Activate()
{
    // don't forget to call base class
    FENLConstraint::Activate();
    m_lc.Activate();
}

//-----------------------------------------------------------------------------
bool FEConstraintImmersedBody::Init()
{
    // initialize surface
//    m_surf.Init();
    
    // we assume that the first domain is the "fluid" domain
    FEMesh& mesh = GetMesh();
    FEDomain& dom = mesh.Domain(0);
    // tag DOFs of all nodes of this domain
    m_nodeBCs.assign(dom.Nodes(),vector<int>(4));
    for (int i=0; i<dom.Nodes(); ++i) {
        int n = dom.NodeIndex(i);
        FENode& node = dom.Node(i);
        for (int k=0; k<3; ++k)
            m_nodeBCs[i][k] = node.get_bc(m_dofW[k]);
        m_nodeBCs[i][3] = node.get_bc(m_dofEF[0]);
    }
    
    // create an FEEdgeList for this domain
    FEEdgeList edgeList;
    edgeList.Create(&dom);
    FEElementEdgeList elemEdgeList;
    elemEdgeList.Create(dom, edgeList);
    
/*
    // evaluate the nodal normals
    int N = m_surf.Nodes(), jp1, jm1;
    vec3d y[FEElement::MAX_NODES], n;
    vector<vec3d>m_nn(N,vec3d(0,0,0));
    
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

    // get the dofs
    int dof_wx = GetDOFIndex("wx");
    int dof_wy = GetDOFIndex("wy");
    int dof_wz = GetDOFIndex("wz");

    // create linear constraints
    // for a frictionless wall the constraint on (vx, vy, vz) is
    // nx*vx + ny*vy + nz*vz = 0
    for (int i=0; i<N; ++i) {
        FEAugLagLinearConstraint* pLC = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        for (int j=0; j<3; ++j) {
            FENode& node = m_surf.Node(i);
            switch (j) {
            case 0: pLC->AddDOF(node.GetID(), dof_wx, m_nn[i].x); break;
            case 1: pLC->AddDOF(node.GetID(), dof_wy, m_nn[i].y); break;
            case 2: pLC->AddDOF(node.GetID(), dof_wz, m_nn[i].z); break;
            default:
                break;
            }
        }
        // add the linear constraint to the system
        m_lc.add(pLC);
    }

    return m_lc.Init();*/
    return     m_surf.Init();
}

//-----------------------------------------------------------------------------
void FEConstraintImmersedBody::Update()
{
    // get the current edge intersection list
    GetIntersectedEdges();
    
    // we assume that the first domain is the "fluid" domain
    FEDomain& dom = GetMesh().Domain(0);
    // for fluid nodes inside the immersed body, prescribe the degrees of freedom to 0
    for (int i=0; i<m_nodetag.size(); ++i) {
        if (m_nodetag[i] == 2) {
            FENode& node = dom.Node(i);
            for (int k=0; k<3; ++k) {
                if (m_nodeBCs[i][k] == DOF_OPEN) {
                    node.set_bc(m_dofW[k], DOF_PRESCRIBED);
                    node.set(m_dofW[k], 0);
                }
            }
//            if (m_nodeBCs[i][3] == DOF_OPEN) {
//                node.set_bc(m_dofEF[0], DOF_PRESCRIBED);
//                node.set(m_dofEF[0], 0);
//            }
        }
/*        else if (m_nodetag[i] == 0) {
            FENode& node = dom.Node(i);
            for (int k=0; k<3; ++k) {
                if (m_nodeBCs[i][k] == DOF_OPEN)
                    node.set_bc(m_dofW[k], DOF_OPEN);
            }
        }*/
    }
}

//-----------------------------------------------------------------------------
void FEConstraintImmersedBody::GetIntersectedEdges()
{
    FEMesh& mesh = GetMesh();
    // we assume that the first domain is the "fluid" domain
    FEDomain& dom = mesh.Domain(0);
    m_EL = FindIntersectedEdges(&dom, &m_surf, m_nodetag);
}

//-----------------------------------------------------------------------------
void FEConstraintImmersedBody::Serialize(DumpStream& ar) { m_lc.Serialize(ar); }
void FEConstraintImmersedBody::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FEConstraintImmersedBody::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FEConstraintImmersedBody::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FEConstraintImmersedBody::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
