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
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEConstraintNormalFlow::Activate()
{
    // don't forget to call base class
    FENLConstraint::Activate();
    m_lc.Activate();
}

//-----------------------------------------------------------------------------
bool FEConstraintNormalFlow::Init()
{
    // initialize surface
    m_surf.Init();
    
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
    for (int i=0; i<N; ++i) {

        //  (1-nx^2)*vx - nx*ny*vy - nx*nz*vz = 0
        FEAugLagLinearConstraint* pLC0 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        for (int j=0; j<3; ++j) {
            FENode& node = m_surf.Node(i);
            switch (j) {
            case 0: pLC0->AddDOF(node.GetID(), dof_wx, 1.0 - m_nn[i].x * m_nn[i].x); break;
            case 1: pLC0->AddDOF(node.GetID(), dof_wy,     - m_nn[i].x * m_nn[i].y); break;
            case 2: pLC0->AddDOF(node.GetID(), dof_wz,     - m_nn[i].x * m_nn[i].z); break;
            }
        }
        // add the linear constraint to the system
        m_lc.add(pLC0);
        
        // -nx*ny*vx + (1-ny^2)*vy - ny*nz*vz = 0
        FEAugLagLinearConstraint* pLC1 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        for (int j=0; j<3; ++j) {
            FENode& node = m_surf.Node(i);
            switch (j) {
                case 0: pLC1->AddDOF(node.GetID(), dof_wx,     -m_nn[i].y * m_nn[i].x); break;
                case 1: pLC1->AddDOF(node.GetID(), dof_wy, 1.0 -m_nn[i].y * m_nn[i].y); break;
                case 2: pLC1->AddDOF(node.GetID(), dof_wz,     -m_nn[i].y * m_nn[i].z); break;
                default:
                    break;
            }
        }
        // add the linear constraint to the system
        m_lc.add(pLC1);
        
        // -nx*nz*vx - ny*nz*vy + (1-nz^2)*vz = 0
        FEAugLagLinearConstraint* pLC2 = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
        for (int j=0; j<3; ++j) {
            FENode node = m_surf.Node(i);
            switch (j) {
                case 0: pLC2->AddDOF(node.GetID(), dof_wx,    -m_nn[i].z * m_nn[i].x); break;
                case 1: pLC2->AddDOF(node.GetID(), dof_wy,    -m_nn[i].z * m_nn[i].y); break;
                case 2: pLC2->AddDOF(node.GetID(), dof_wz, 1.0-m_nn[i].z * m_nn[i].z); break;
                default:
                    break;
            }
        }
        // add the linear constraint to the system
        m_lc.add(pLC2);
    }
    
    return true;
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
