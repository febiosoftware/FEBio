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
#include "FEConstraintFrictionlessWall.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FESurface.h>

BEGIN_FECORE_CLASS(FEConstraintFrictionlessWall, FESurfaceConstraint)
    ADD_PARAMETER(m_lc.m_laugon, "laugon");
    ADD_PARAMETER(m_lc.m_tol, "tol");
    ADD_PARAMETER(m_lc.m_eps, "penalty");
    ADD_PARAMETER(m_lc.m_rhs, "rhs");
    ADD_PARAMETER(m_lc.m_naugmin, "minaug");
    ADD_PARAMETER(m_lc.m_naugmax, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEConstraintFrictionlessWall::FEConstraintFrictionlessWall(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem), m_lc(pfem)
{
    m_binit = false;
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEConstraintFrictionlessWall::Activate()
{
    // don't forget to call base class
    FESurfaceConstraint::Activate();
    
    if (m_binit == false) {
        
        // evaluate the nodal normals
        m_surf.UpdateNodeNormals();
        
        // get the dofs
        int dof_wx = GetDOFIndex("wx");
        int dof_wy = GetDOFIndex("wy");
        int dof_wz = GetDOFIndex("wz");

        // create linear constraints
        // for a frictionless wall the constraint on (vx, vy, vz) is
        // nx*vx + ny*vy + nz*vz = 0
        for (int i=0; i<m_surf.Nodes(); ++i) {
            FEAugLagLinearConstraint* pLC = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
            FENode& node = m_surf.Node(i);
            vec3d nn = m_surf.NodeNormal(i);
            pLC->AddDOF(node.GetID(), dof_wx, nn.x);
            pLC->AddDOF(node.GetID(), dof_wy, nn.y);
            pLC->AddDOF(node.GetID(), dof_wz, nn.z);
            // add the linear constraint to the system
            m_lc.add(pLC);
        }
        m_lc.Init();
        m_lc.Activate();
        m_binit = true;
    }
}

//-----------------------------------------------------------------------------
bool FEConstraintFrictionlessWall::Init()
{
    // initialize surface
    return m_surf.Init();
}

//-----------------------------------------------------------------------------
void FEConstraintFrictionlessWall::Serialize(DumpStream& ar) { m_lc.Serialize(ar); }
void FEConstraintFrictionlessWall::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FEConstraintFrictionlessWall::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FEConstraintFrictionlessWall::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FEConstraintFrictionlessWall::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
int FEConstraintFrictionlessWall::InitEquations(int neq) { return m_lc.InitEquations(neq); }
void FEConstraintFrictionlessWall::Update(const std::vector<double>& Ui, const std::vector<double>& ui) { m_lc.Update(Ui, ui); }
void FEConstraintFrictionlessWall::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui) { m_lc.UpdateIncrements(Ui, ui); }
void FEConstraintFrictionlessWall::PrepStep() { m_lc.PrepStep(); }
