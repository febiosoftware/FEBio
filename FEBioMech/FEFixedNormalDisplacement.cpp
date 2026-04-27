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
#include "FEFixedNormalDisplacement.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFixedNormalDisplacement, FESurfaceConstraint)
    ADD_PARAMETER(m_lc.m_laugon, "laugon");
    ADD_PARAMETER(m_lc.m_tol, "tol");
    ADD_PARAMETER(m_lc.m_eps, "penalty");
    ADD_PARAMETER(m_lc.m_rhs, "rhs");
    ADD_PARAMETER(m_lc.m_naugmin, "minaug");
    ADD_PARAMETER(m_lc.m_naugmax, "maxaug");

    ADD_PARAMETER(m_bshellb , "shell_bottom");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFixedNormalDisplacement::FEFixedNormalDisplacement(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem), m_lc(pfem)
{
    m_binit = false;
    m_bshellb = false;
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEFixedNormalDisplacement::Activate()
{
    // don't forget to call base class
    FESurfaceConstraint::Activate();

    if (m_binit == false)
    {
        FEModel& fem = *GetFEModel();
        DOFS& dofs = fem.GetDOFS();

        // evaluate the nodal normals
        m_surf.UpdateNodeNormals();

        // create linear constraints
        // for a symmetry plane the constraint on (ux, uy, uz) is
        // nx*ux + ny*uy + nz*uz = 0
        if (m_bshellb == false) {
            for (int i = 0; i < m_surf.Nodes(); ++i) {
                FENode& node = m_surf.Node(i);
                if ((node.HasFlags(FENode::EXCLUDE) == false) && (node.m_rid == -1)) {
					int nid = m_surf.NodeIndex(i) + 1; // NOTE: the linear constraint set still assumes 1-based indexing!
                    vec3d nn = m_surf.NodeNormal(i);
                    FEAugLagLinearConstraint* pLC = fecore_alloc(FEAugLagLinearConstraint, &fem);
                    pLC->AddDOF(nid, dofs.GetDOF("x"), nn.x);
                    pLC->AddDOF(nid, dofs.GetDOF("y"), nn.y);
                    pLC->AddDOF(nid, dofs.GetDOF("z"), nn.z);
                    // add the linear constraint to the system
                    m_lc.add(pLC);
                }
            }
        }
        else {
            for (int i = 0; i < m_surf.Nodes(); ++i) {
                FENode& node = m_surf.Node(i);
                if ((node.HasFlags(FENode::EXCLUDE) == false) && (node.m_rid == -1)) {
                    vec3d nn = m_surf.NodeNormal(i);
					int nid = m_surf.NodeIndex(i) + 1; // NOTE: the linear constraint set still assumes 1-based indexing!
					FEAugLagLinearConstraint* pLC = fecore_alloc(FEAugLagLinearConstraint, &fem);
                    pLC->AddDOF(nid, dofs.GetDOF("sx"), nn.x);
                    pLC->AddDOF(nid, dofs.GetDOF("sy"), nn.y);
                    pLC->AddDOF(nid, dofs.GetDOF("sz"), nn.z);
                    // add the linear constraint to the system
                    m_lc.add(pLC);
                }
            }
        }

        m_lc.Init();
        m_lc.Activate();
        m_binit = true;
    }
}

//-----------------------------------------------------------------------------
bool FEFixedNormalDisplacement::Init()
{
    // initialize surface
    return m_surf.Init();
}

//-----------------------------------------------------------------------------
void FEFixedNormalDisplacement::Serialize(DumpStream& ar) { m_lc.Serialize(ar); }
void FEFixedNormalDisplacement::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FEFixedNormalDisplacement::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FEFixedNormalDisplacement::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FEFixedNormalDisplacement::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
int FEFixedNormalDisplacement::InitEquations(int neq) { return m_lc.InitEquations(neq); }
void FEFixedNormalDisplacement::Update(const std::vector<double>& Ui, const std::vector<double>& ui) { m_lc.Update(Ui, ui); }
void FEFixedNormalDisplacement::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui) { m_lc.UpdateIncrements(Ui, ui); }
void FEFixedNormalDisplacement::PrepStep() { m_lc.PrepStep(); }
