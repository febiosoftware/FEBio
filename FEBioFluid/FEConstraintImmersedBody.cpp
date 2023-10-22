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
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>

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
    m_breset = true;
    static int ns = 0;
    // TODO: Can this be done in Init, since  there is no error checking
    if (pfem)
    {
        m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
        m_dofEF.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION));
        m_nodalLoad = fecore_alloc(FENodalDOFLoad, pfem);
        m_nodalLoad->SetDOF(m_dofEF[0]);
        m_nodalLoad->SetLoad(0);
        char name[256];
        snprintf(name,256,"ImmersedBodyZeroFluidFlux%d",ns++);
        m_nodalLoad->SetName(name);
        FENodeSet* nset = new FENodeSet(GetFEModel());
        m_nodalLoad->SetNodeSet(nset);
        pfem->AddModelLoad(m_nodalLoad);

        FEMesh& mesh = GetMesh();
        int nbc = pfem->BoundaryConditions();
        for (int i=0; i<nbc; ++i)
        {
            FEBoundaryCondition* bc = pfem->BoundaryCondition(i);
            FEPrescribedNodeSet* nset = dynamic_cast<FEPrescribedNodeSet*>(bc);
            if (nset) {
                if (nset->GetName() == string("PrescribedFluidVelocityX")) m_pbcwx = nset;
                else if (nset->GetName() == string("PrescribedFluidVelocityY")) m_pbcwy = nset;
                else if (nset->GetName() == string("PrescribedFluidVelocityZ")) m_pbcwz = nset;
                else if (nset->GetName() == string("PrescribedFluidDilatation")) m_pbcef = nset;
            }
        }
        if (m_pbcwx) m_pbcwx->GetNodeSet()->Clear();
        if (m_pbcwy) m_pbcwy->GetNodeSet()->Clear();
        if (m_pbcwz) m_pbcwz->GetNodeSet()->Clear();
        if (m_pbcef) m_pbcef->GetNodeSet()->Clear();
    }
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEConstraintImmersedBody::Activate()
{
    // don't forget to call base class
    FENLConstraint::Activate();
    m_lc.Activate();
    m_nodalLoad->Activate();
}

//-----------------------------------------------------------------------------
bool FEConstraintImmersedBody::Init()
{
/*
    // create an FEEdgeList for this domain
    FEEdgeList edgeList;
    edgeList.Create(&dom);
    FEElementEdgeList elemEdgeList;
    elemEdgeList.Create(dom, edgeList);
*/
    return     m_surf.Init();
}

//-----------------------------------------------------------------------------
void FEConstraintImmersedBody::PrepStep()
{
    // get the current edge intersection list
    GetIntersectedEdges();
    
    // we assume that the first domain is the "fluid" domain
    FEMesh& mesh = GetMesh();
    FEDomain& dom = mesh.Domain(0);
/*
    if (m_breset) {
        // tag DOFs of all nodes of this domain
        m_nodeBCs.assign(dom.Nodes(),vector<int>(4));
        for (int i=0; i<dom.Nodes(); ++i) {
            FENode& node = dom.Node(i);
            for (int k=0; k<3; ++k)
                m_nodeBCs[i][k] = node.get_bc(m_dofW[k]);
            m_nodeBCs[i][3] = node.get_bc(m_dofEF[0]);
        }
        m_breset = false;
    }
*/
    // for fluid nodes inside the immersed body, prescribe the degrees of freedom
    if (m_pbcwx) m_pbcwx->GetNodeSet()->Clear();
    if (m_pbcwy) m_pbcwy->GetNodeSet()->Clear();
    if (m_pbcwz) m_pbcwz->GetNodeSet()->Clear();
    if (m_pbcef) m_pbcef->GetNodeSet()->Clear();
    for (int i=0; i<m_nodetag.size(); ++i) {
//        if (m_nodetag[i] > 0) {
        if (m_nodetag[i] == 2) {
            int nid = dom.NodeIndex(i);
            // for now, set the velocity DOFs to zero
            if (m_pbcwx) m_pbcwx->GetNodeSet()->Add(nid);
            if (m_pbcwy) m_pbcwy->GetNodeSet()->Add(nid);
            if (m_pbcwz) m_pbcwz->GetNodeSet()->Add(nid);
            if (m_pbcef) m_pbcef->GetNodeSet()->Add(nid);
        }
    }
    if (m_pbcwx) { m_pbcwx->Activate(); m_pbcwx->Repair(); }
    if (m_pbcwy) { m_pbcwy->Activate(); m_pbcwy->Repair(); }
    if (m_pbcwz) { m_pbcwz->Activate(); m_pbcwz->Repair(); }
    if (m_pbcef) { m_pbcef->Activate(); m_pbcef->Repair(); }

    // prescribe linear constraints on intersected edges
    m_lc.m_LC.clear();
    for (int i=0; i<m_EL.Edges(); ++i) {
        FEEdgeList::EDGE EL = m_EL.Edge(i);
        FENode& node0 = dom.Node(EL.node[0]);
        FENode& node1 = dom.Node(EL.node[1]);
        int nid0 = node0.GetID();
        int nid1 = node1.GetID();
        double g = EL.tag;

        // create linear constraints
        FEAugLagLinearConstraint* pLC0 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        pLC0->AddDOF(nid0, m_dofW[0], 1-g);
        pLC0->AddDOF(nid1, m_dofW[0],   g);
        pLC0->SetRHS(0.0);
        FEAugLagLinearConstraint* pLC1 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        pLC1->AddDOF(nid0, m_dofW[1], 1-g);
        pLC1->AddDOF(nid1, m_dofW[1],   g);
        pLC1->SetRHS(0.0);
        FEAugLagLinearConstraint* pLC2 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        pLC2->AddDOF(nid0, m_dofW[2], 1-g);
        pLC2->AddDOF(nid1, m_dofW[2],   g);
        pLC2->SetRHS(0.0);
        // add the linear constraint to the system
        m_lc.add(pLC0);
        m_lc.add(pLC1);
        m_lc.add(pLC2);
    }

    // populate the node set for the nodal load
    FENodeSet* nset = m_nodalLoad->GetNodeSet();
    nset->Clear();
    for (int i=0; i<m_nodetag[i]; ++i) {
        if ((m_nodetag[i] == 1) || (m_nodetag[i] == -1)) {
            int nid = dom.NodeIndex(i);
            nset->Add(nid);
        }
    }

    // force mesh update
    GetFEModel()->SetMeshUpdateFlag(true);
}

//-----------------------------------------------------------------------------
void FEConstraintImmersedBody::GetIntersectedEdges()
{
    FEMesh& mesh = GetMesh();
    // we assume that the first domain is the "fluid" domain
    FEDomain& dom = mesh.Domain(0);
    m_EL = FindIntersectedEdges(&dom, &m_surf, m_nodetag, m_edgetag);
}

//-----------------------------------------------------------------------------
void FEConstraintImmersedBody::Serialize(DumpStream& ar) { m_lc.Serialize(ar); }
void FEConstraintImmersedBody::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FEConstraintImmersedBody::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FEConstraintImmersedBody::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FEConstraintImmersedBody::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
