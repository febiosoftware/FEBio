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
    ADD_PARAMETER(m_lc.m_laugon, "laugon");
    ADD_PARAMETER(m_lc.m_tol, "tol");
    ADD_PARAMETER(m_lc.m_eps, "penalty");
    ADD_PARAMETER(m_lc.m_rhs, "rhs");
    ADD_PARAMETER(m_lc.m_naugmin, "minaug");
    ADD_PARAMETER(m_lc.m_naugmax, "maxaug");
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
    
    // get the current edge intersection list
    GetIntersectedEdges();
    
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
    /*
    // evaluate the current average value of the dilatation for nodes with tag = +1
    double ef = 0;
    int in = 0;
    for (int i=0; i<m_nodetag.size(); ++i) {
        if (m_nodetag[i] == 1) {
            FENode& node = dom.Node(i);
            ef += node.get(m_dofEF[0]);
            ++in;
        }
    }
    if (in) ef /= in;*/
    // for fluid nodes inside the immersed body, prescribe the degrees of freedom
    for (int i=0; i<m_nodetag.size(); ++i) {
        if (m_nodetag[i] == 2) {
            FENode& node = dom.Node(i);
            // for now, set the velocity DOFs to zero
            for (int k=0; k<3; ++k) {
                if (m_nodeBCs[i][k] == DOF_OPEN) {
                    node.set_bc(m_dofW[k], DOF_PRESCRIBED);
                    node.set(m_dofW[k], 0);
                }
            }
        }
    }
    
    // prescribe linear constraints on intersected edges
    for (int i=0; i<m_EL.Edges(); ++i) {
        FEEdgeList::EDGE EL = m_EL.Edge(i);
        FENode& node0 = mesh.Node(EL.node[0]);
        FENode& node1 = mesh.Node(EL.node[1]);
        double g = EL.tag;
        
        // create linear constraints
        FEAugLagLinearConstraint* pLC0 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        pLC0->AddDOF(node0.GetID(), m_dofW[0], 1-g);
        pLC0->AddDOF(node1.GetID(), m_dofW[0],   g);
        FEAugLagLinearConstraint* pLC1 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        pLC1->AddDOF(node0.GetID(), m_dofW[1], 1-g);
        pLC1->AddDOF(node1.GetID(), m_dofW[1],   g);
        FEAugLagLinearConstraint* pLC2 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        pLC2->AddDOF(node0.GetID(), m_dofW[2], 1-g);
        pLC2->AddDOF(node1.GetID(), m_dofW[2],   g);
        // add the linear constraint to the system
        m_lc.add(pLC0);
        m_lc.add(pLC1);
        m_lc.add(pLC2);
    }

    /*
    // create an FEEdgeList for this domain
    FEEdgeList edgeList;
    edgeList.Create(&dom);
    FEElementEdgeList elemEdgeList;
    elemEdgeList.Create(dom, edgeList);*/
    
    return     m_surf.Init();
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
