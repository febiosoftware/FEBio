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
#include "FEConstraintImmersedFSIBody.h"
#include "FEBioFluid.h"
#include <FECore/FEElementList.h>
#include <FECore/FEEdgeList.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>

BEGIN_FECORE_CLASS(FEConstraintImmersedFSIBody, FESurfaceConstraint)
    ADD_PARAMETER(m_lc.m_laugon, "laugon");
    ADD_PARAMETER(m_lc.m_tol, "tol");
    ADD_PARAMETER(m_lc.m_eps, "penalty");
    ADD_PARAMETER(m_lc.m_rhs, "rhs");
    ADD_PARAMETER(m_lc.m_naugmin, "minaug");
    ADD_PARAMETER(m_lc.m_naugmax, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEConstraintImmersedFSIBody::FEConstraintImmersedFSIBody(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem), m_lc(pfem), m_dofW(pfem), m_dofEF(pfem)
{
    static int ns = 0;
    // TODO: Can this be done in Init, since  there is no error checking
    if (pfem)
    {
        m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
        m_dofEF.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION));
/*        m_nodalLoad = fecore_alloc(FENodalDOFLoad, pfem);
        m_nodalLoad->SetDOF(m_dofEF[0]);
        m_nodalLoad->SetLoad(0);
        char name[256];
        snprintf(name,256,"ImmersedBodyZeroFluidFlux%d",ns++);
        m_nodalLoad->SetName(name);
        FENodeSet* nset = new FENodeSet(GetFEModel());
        m_nodalLoad->SetNodeSet(nset);
        pfem->AddModelLoad(m_nodalLoad);*/
    }
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEConstraintImmersedFSIBody::Activate()
{
    // don't forget to call base class
    FENLConstraint::Activate();
    m_lc.Activate();
//    m_nodalLoad->Activate();
}

//-----------------------------------------------------------------------------
bool FEConstraintImmersedFSIBody::Init()
{
    FEMesh& mesh = GetMesh();
    // find the surface of the FluidFSITraction by name
    FESurface* FSItsrf = mesh.FindSurface("FluidFSITraction2");
    if (FSItsrf == nullptr) return false;
    
    // get the current edge intersection list
//    GetIntersectedEdges(FSItsrf);
    
    // we assume that the first domain is the "fluid" domain
    FEDomain& fluiddom = mesh.Domain(0);
    // tag DOFs of all nodes of this domain
    m_nodeBCs.assign(fluiddom.Nodes(),vector<int>(4));
    for (int i=0; i<fluiddom.Nodes(); ++i) {
        FENode& node = fluiddom.Node(i);
        for (int k=0; k<3; ++k)
            m_nodeBCs[i][k] = node.get_bc(m_dofW[k]);
        m_nodeBCs[i][3] = node.get_bc(m_dofEF[0]);
    }
    
    // we assume that the second domain is the "solid" domain
    m_nodetag.assign(fluiddom.Nodes(),0);
    FEDomain& soliddom = mesh.Domain(1);
    FESolidDomain* sdom = dynamic_cast<FESolidDomain*>(&soliddom);
    if (sdom == nullptr) return false;
    
    bool reset = true;
    int count = 0;
    for (int i=0; i<fluiddom.Nodes(); ++i) {
        int nidx = fluiddom.NodeIndex(i);
        if (FSItsrf->IsInsideSurface(nidx, reset,0.001)) {
            m_nodetag[i] = 2;
            ++count;
        }
    }

    // for each node on the FSI surface check which solid element
    // of domain dom encloses it
    sdom = dynamic_cast<FESolidDomain*>(&fluiddom);
    if (sdom == nullptr) return false;
    vector<FESolidElement*> sel;
    vector<vec3d> selrs;
    for (int i=0; i<m_surf.Nodes(); ++i) {
        FENode& node = m_surf.Node(i);
        vec3d p = node.m_rt;
        double r[3];
        FESolidElement* s = sdom->FindElement(p, r);
        if (s) {
            sel.push_back(s);
            selrs.push_back(vec3d(r[0],r[1],r[2]));
        }
    }
    if (sel.size() != m_surf.Nodes())
        return false;

    // for fluid nodes inside the immersed solid body, prescribe the degrees of freedom
    for (int i=0; i<m_nodetag.size(); ++i) {
//        if (m_nodetag[i] == 2) {
        if (m_nodetag[i] > 0) {
            FENode& node = fluiddom.Node(i);
            // for now, set the velocity DOFs to zero
            for (int k=0; k<3; ++k) {
                if (m_nodeBCs[i][k] == DOF_OPEN) {
                    node.set_bc(m_dofW[k], DOF_PRESCRIBED);
                    node.set(m_dofW[k], 0);
                }
            }
            // prevent dilatation growing out of bounds inside solid domain
            if (m_nodeBCs[i][3] == DOF_OPEN) {
                node.set_bc(m_dofEF[0], DOF_PRESCRIBED);
                node.set(m_dofEF[0], node.get_prev(m_dofEF[0]));
            }
        }
    }

    // prescribe linear constraints on elements
    double kappa = 1e2;
    for (int i=0; i<m_surf.Nodes(); ++i) {
        FENode& node = m_surf.Node(i);
        FESolidElement* se = sel[i];
        vec3d rs = selrs[i];
        double H[FEElement::MAX_NODES];
        se->shape_fnc(H, rs.x, rs.y, rs.z);
        // create linear constraints
        FEAugLagLinearConstraint* pLCX = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        FEAugLagLinearConstraint* pLCY = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        FEAugLagLinearConstraint* pLCZ = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        FEAugLagLinearConstraint* pLCE = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        pLCX->AddDOF(node.GetID(), m_dofW[0], -1);
        pLCY->AddDOF(node.GetID(), m_dofW[1], -1);
        pLCZ->AddDOF(node.GetID(), m_dofW[2], -1);
        pLCE->AddDOF(node.GetID(), m_dofEF[0],-1*kappa);
        for (int j=0; j<se->Nodes(); ++j) {
            FENode& enode = mesh.Node(se->m_node[j]);
            pLCX->AddDOF(enode.GetID(), m_dofW[0], H[j]);
            pLCY->AddDOF(enode.GetID(), m_dofW[1], H[j]);
            pLCZ->AddDOF(enode.GetID(), m_dofW[2], H[j]);
            pLCE->AddDOF(enode.GetID(), m_dofEF[0], H[j]*kappa);
        }
        // add the linear constraint to the system
        m_lc.add(pLCX);
        m_lc.add(pLCY);
        m_lc.add(pLCZ);
        m_lc.add(pLCE);
    }
/*
    // populate the node set for the nodal load
    FENodeSet* nset = m_nodalLoad->GetNodeSet();
    for (int i=0; i<m_nodetag[i]; ++i) {
        if ((m_nodetag[i] == 1) || (m_nodetag[i] == -1)) {
            int nid = fluiddom.NodeIndex(i);
            nset->Add(nid);
        }
    }
*/
    return m_surf.Init();
}

//-----------------------------------------------------------------------------
void FEConstraintImmersedFSIBody::GetIntersectedEdges(FESurface* surf)
{
    FEMesh& mesh = GetMesh();
    // we assume that the first domain is the "fluid" domain
    FEDomain& dom = mesh.Domain(0);
    m_EL = FindIntersectedEdges(&dom, surf, m_nodetag, m_edgetag);
}

//-----------------------------------------------------------------------------
void FEConstraintImmersedFSIBody::Serialize(DumpStream& ar) { m_lc.Serialize(ar); }
void FEConstraintImmersedFSIBody::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FEConstraintImmersedFSIBody::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FEConstraintImmersedFSIBody::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FEConstraintImmersedFSIBody::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
