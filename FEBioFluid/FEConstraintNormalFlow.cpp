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
#include "FEConstraintNormalFlow.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEMesh.h>

BEGIN_FECORE_CLASS(FEConstraintNormalFlow, FESurfaceConstraint)
    ADD_PARAMETER(m_lc.m_laugon, "laugon");
    ADD_PARAMETER(m_lc.m_tol, "tol");
    ADD_PARAMETER(m_lc.m_eps, "penalty");
    ADD_PARAMETER(m_lc.m_rhs, "rhs");
    ADD_PARAMETER(m_lc.m_naugmin, "minaug");
    ADD_PARAMETER(m_lc.m_naugmax, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEConstraintNormalFlow::FEConstraintNormalFlow(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem), m_lc(pfem)
{
    m_binit = false;
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEConstraintNormalFlow::Activate()
{
    // don't forget to call base class
    FESurfaceConstraint::Activate();
    
    if (m_binit == false)
    {
        m_surf.UpdateNodeNormals();
        
        // get the dof indices
        int dof_wx = GetDOFIndex("wx");
        int dof_wy = GetDOFIndex("wy");
        int dof_wz = GetDOFIndex("wz");
        
        // create linear constraints
        // for a surface with zero tangential velocity the constraints
        // on (vx, vy, vz) are
        //  (1-nx^2)*vx - nx*ny*vy - nx*nz*vz = 0
        // -nx*ny*vx + (1-ny^2)*vy - ny*nz*vz = 0
        // -nx*nz*vx - ny*nz*vy + (1-nz^2)*vz = 0
        for (int i=0; i<m_surf.Nodes(); ++i) {
            FENode& node = m_surf.Node(i);
            if ((node.HasFlags(FENode::EXCLUDE) == false) && (node.m_rid == -1)) {
                vec3d nn = m_surf.NodeNormal(i);
                
                //  (1-nx^2)*vx - nx*ny*vy - nx*nz*vz = 0
                FEAugLagLinearConstraint* pLC0 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                pLC0->AddDOF(node.GetID(), dof_wx, 1.0 - nn.x * nn.x);
                pLC0->AddDOF(node.GetID(), dof_wy,     - nn.x * nn.y);
                pLC0->AddDOF(node.GetID(), dof_wz,     - nn.x * nn.z);
                // add the linear constraint to the system
                if ((pLC0->m_dof[0]->m_val != 0) || (pLC0->m_dof[1]->m_val != 0) || (pLC0->m_dof[2]->m_val != 0))
                    m_lc.add(pLC0);
                
                // -nx*ny*vx + (1-ny^2)*vy - ny*nz*vz = 0
                FEAugLagLinearConstraint* pLC1 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                pLC1->AddDOF(node.GetID(), dof_wx,     -nn.y * nn.x);
                pLC1->AddDOF(node.GetID(), dof_wy, 1.0 -nn.y * nn.y);
                pLC1->AddDOF(node.GetID(), dof_wz,     -nn.y * nn.z);
                // add the linear constraint to the system
                if ((pLC1->m_dof[0]->m_val != 0) || (pLC1->m_dof[1]->m_val != 0) || (pLC1->m_dof[2]->m_val != 0))
                    m_lc.add(pLC1);
                
                // -nx*nz*vx - ny*nz*vy + (1-nz^2)*vz = 0
                FEAugLagLinearConstraint* pLC2 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                pLC2->AddDOF(node.GetID(), dof_wx,    -nn.z * nn.x);
                pLC2->AddDOF(node.GetID(), dof_wy,    -nn.z * nn.y);
                pLC2->AddDOF(node.GetID(), dof_wz, 1.0-nn.z * nn.z);
                // add the linear constraint to the system
                if ((pLC2->m_dof[0]->m_val != 0) || (pLC2->m_dof[1]->m_val != 0) || (pLC2->m_dof[2]->m_val != 0))
                    m_lc.add(pLC2);
            }
        }
        m_lc.Init();
        m_lc.Activate();
        m_binit = true;
    }
}

//-----------------------------------------------------------------------------
bool FEConstraintNormalFlow::Init()
{
    // initialize surface
    return m_surf.Init();
}

//-----------------------------------------------------------------------------
void FEConstraintNormalFlow::Serialize(DumpStream& ar)
{
	FESurfaceConstraint::Serialize(ar);
    m_lc.Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_surf;
}

//-----------------------------------------------------------------------------
void FEConstraintNormalFlow::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FEConstraintNormalFlow::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FEConstraintNormalFlow::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FEConstraintNormalFlow::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
int FEConstraintNormalFlow::InitEquations(int neq) { return m_lc.InitEquations(neq); }
void FEConstraintNormalFlow::Update(const std::vector<double>& Ui, const std::vector<double>& ui) { m_lc.Update(Ui, ui); }
void FEConstraintNormalFlow::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui) { m_lc.UpdateIncrements(Ui, ui); }
void FEConstraintNormalFlow::PrepStep() { m_lc.PrepStep(); }
