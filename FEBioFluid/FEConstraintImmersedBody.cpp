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
#include "FEBioFSI.h"
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
FEConstraintImmersedBody::FEConstraintImmersedBody(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem), m_lc(pfem), m_dofW(pfem), m_dofEF(pfem), m_dofV(pfem)
{
    m_breset = true;
    static int ns = 0;
    // TODO: Can this be done in Init, since  there is no error checking
    if (pfem)
    {
        m_dofW.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_VELOCITY));
        m_dofEF.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION));
        m_dofV.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::VELOCITY));
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
    // we assume that the first domain is the "fluid" domain
    // determine which elements connect to each node in this domain
    FEMesh& mesh = GetMesh();
    FEDomain& dom = mesh.Domain(0);
    m_nodal_elems.resize(dom.Nodes());
    for (int i=0; i<dom.Elements(); ++i) {
        FEElement& el = dom.ElementRef(i);
        for (int j=0; j<el.Nodes(); ++j)
            m_nodal_elems[el.m_lnode[j]].push_back(i);
    }
    
    return     m_surf.Init();
}

//-----------------------------------------------------------------------------
void FEConstraintImmersedBody::PrepStep()
{
    // get the current edge intersection list
    GetIntersectedEdges();
    IntersectionSurface();
    
    // we assume that the first domain is the "fluid" domain
    FEMesh& mesh = GetMesh();
    FEDomain& dom = mesh.Domain(0);

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
    if (m_pbcwx) { m_pbcwx->Activate(); m_pbcwx->Repair(); m_pbcwx->SetRelativeFlag(true); }
    if (m_pbcwy) { m_pbcwy->Activate(); m_pbcwy->Repair(); m_pbcwy->SetRelativeFlag(true); }
    if (m_pbcwz) { m_pbcwz->Activate(); m_pbcwz->Repair(); m_pbcwz->SetRelativeFlag(true); }
    if (m_pbcef) { m_pbcef->Activate(); m_pbcef->Repair(); m_pbcef->SetRelativeFlag(true); }

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
        pLC0->SetRHS(-m_vs[i].x);
        FEAugLagLinearConstraint* pLC1 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        pLC1->AddDOF(nid0, m_dofW[1], 1-g);
        pLC1->AddDOF(nid1, m_dofW[1],   g);
        pLC1->SetRHS(-m_vs[i].y);
        FEAugLagLinearConstraint* pLC2 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        pLC2->AddDOF(nid0, m_dofW[2], 1-g);
        pLC2->AddDOF(nid1, m_dofW[2],   g);
        pLC2->SetRHS(-m_vs[i].z);
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
FESurface* FEConstraintImmersedBody::IntersectionSurface()
{
    // we assume that the first domain is the "fluid" domain
    FEMesh& mesh = GetMesh();
    FEDomain& dom = mesh.Domain(0);
    FETimeInfo& tp = GetFEModel()->GetTime();
    
    // reset solid velocity vectors
    int NE = m_EL.Edges();
    m_vs.assign(NE,vec3d(0,0,0));
    m_vsn.assign(NE,0.0);
    
    // surface normal at the intersection points
    vec3d normal(0,0,0);
    // fluid elements thad include an intersected edge
    vector< vector<int>> edge_elems(m_EL.Edges());
    
    for (int i=0; i<m_EL.Edges(); ++i)
    {
        FEEdgeList::EDGE edge = m_EL.Edge(i);
        // Find the normal to the solid face at the edge intersection point
        FESurfaceElement& sel = m_surf.Element(edge.selid);
        normal = m_surf.SurfaceNormal(sel, edge.rs.x(), edge.rs.y());
        // find the velocity of the solid face at the edge intersection point
        vec3d vs[FEElement::MAX_NODES];
        for (int j=0; j<sel.Nodes(); ++j) {
            FENode& node = m_surf.Node(sel.m_lnode[j]);
            vs[j] = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2])*tp.alphaf + node.get_vec3d_prev(m_dofV[0], m_dofV[1], m_dofV[2])*(1-tp.alphaf);
        }
        m_vs[i] = sel.eval(vs, edge.rs.x(), edge.rs.y());
        m_vsn[i] = m_vs[i]*normal;
        
        // Find the elements that contain this edge
        int n1 = edge.node[0];
        int n2 = edge.node[1];
        for (int j=0; j<m_nodal_elems[n1].size(); ++j) {
            for (int k=0; k<m_nodal_elems[n2].size(); ++k) {
                if (m_nodal_elems[n1][j] == m_nodal_elems[n2][k])
                    edge_elems[i].push_back(m_nodal_elems[n1][j]);
            }
        }
    }
    
    // Generate a list of elements that contain the integration points of the immersed surface, and their parametric coordinates
    m_int_elem.resize(0);
    m_int_rst.resize(0);
    m_int_selid.resize(0);
    m_int_semp.resize(0);
    for (int i=0; i<m_EL.Edges(); ++i)
    {
        FEEdgeList::EDGE edge = m_EL.Edge(i);
        // Find the normal to the solid face at the edge intersection point
        FESurfaceElement& sel = m_surf.Element(edge.selid);
        for (int j=0; j<sel.GaussPoints(); ++j) {
            FEMaterialPoint& pt = *(sel.GetMaterialPoint(j));
            vec3d x = pt.m_rt;
            for (int k=0; k<edge_elems[i].size(); ++k) {
                FEElement* e = &dom.ElementRef(edge_elems[i][k]);
                FESolidElement* el = dynamic_cast<FESolidElement*>(e);
                if (el) {
                    FEBoundingBox box;
                    for (int l=0; l<el->Nodes(); ++l) {
                        FENode& node = dom.Node(el->m_lnode[l]);
                        box.add(node.m_rt);
                    }
                    if (box.IsInside(x)) {
                        double r = 0, s = 0, t = 0;
                        if (FindPoint(*el, x, r, s, t)) {
                            m_int_elem.push_back(el);
                            vec3d rst(r,s,t);
                            m_int_rst.push_back(rst);
                            m_int_selid.push_back(edge.selid);
                            m_int_semp.push_back(j);
                        }
                    }
                }
            }
        }
    }
}


//-----------------------------------------------------------------------------
bool FEConstraintImmersedBody::FindPoint(FESolidElement& el, const vec3d x, double& r, double& s, double& t)
{
    // we assume that the first domain is the "fluid" domain
    FEMesh& mesh = GetMesh();
    FEDomain& dom = mesh.Domain(0);
    
    bool    found = false;
    double  errrel = 1e-6;
    double  errabs = 1e-9;
    int     maxiter = 20;
    double  tol = 1e-5;
    
    bool convgd = false;
    int iter = 0;
    int neln = el.Nodes();
    double *N, *Nr, *Ns, *Nt;
    vec3d f(0,0,0);
    vec3d df[3];
    mat3d gradf;
    do {
        el.shape_fnc(N,r,s,t);
        el.shape_deriv(Nr, Ns, Nt, r, s, t);
        f = -x;
        df[0] = df[1] = df[2] = vec3d(0,0,0);
        for (int i=0; i<neln; ++i) {
            FENode& node = dom.Node(el.m_lnode[i]);
            f += node.m_rt*N[i];
            df[0] += node.m_rt*Nr[i];
            df[1] += node.m_rt*Ns[i];
            df[2] += node.m_rt*Nt[i];
        }
        mat3d gradf(df[0].x, df[1].x, df[2].x,
                    df[0].y, df[1].y, df[2].y,
                    df[0].z, df[1].z, df[2].z);
        vec3d dr = -(gradf.inverse()*f);
        r += dr.x; s += dr.y; t += dr.z;
        double rmag = sqrt(r*r+s*s+t*t);
        double drmag = dr.Length();
        double fmag = f.Length();
        if (drmag <= errrel*rmag) convgd = true;
        if (fmag <= errabs) convgd = true;
        ++iter;
    } while (!convgd && (iter < maxiter));
    
    if (convgd) {
        if ((r >= -1-tol) && (r <= 1+tol)
            && (s >= -1-tol) && (s <= 1+tol)
            && (t >= -1-tol) && (t <= 1+tol))
            found = true;
    }
    else found = false;
    
    return found;
}


//-----------------------------------------------------------------------------
void FEConstraintImmersedBody::Serialize(DumpStream& ar) { m_lc.Serialize(ar); }
void FEConstraintImmersedBody::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FEConstraintImmersedBody::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FEConstraintImmersedBody::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FEConstraintImmersedBody::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
