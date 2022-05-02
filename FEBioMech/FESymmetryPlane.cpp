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
#include "FESymmetryPlane.h"
#include "FECore/FECoreKernel.h"
//#include <FECore/FEModel.h>

BEGIN_FECORE_CLASS(FESymmetryPlane, FESurfaceConstraint)
    ADD_PARAMETER(m_lc.m_laugon, "laugon");
    ADD_PARAMETER(m_lc.m_tol, "tol");
    ADD_PARAMETER(m_lc.m_eps, "penalty");
    ADD_PARAMETER(m_lc.m_rhs, "rhs");
    ADD_PARAMETER(m_lc.m_naugmin, "minaug");
    ADD_PARAMETER(m_lc.m_naugmax, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESymmetryPlane::FESymmetryPlane(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem), m_lc(pfem)
{
	m_binit = false;
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FESymmetryPlane::Activate()
{
    // don't forget to call base class
    FESurfaceConstraint::Activate();

    if (m_binit == false)
    {
        // get the dof indices
        int dofX = GetDOFIndex("x");
        int dofY = GetDOFIndex("y");
        int dofZ = GetDOFIndex("z");
        int dofSX = GetDOFIndex("sx");
        int dofSY = GetDOFIndex("sy");
        int dofSZ = GetDOFIndex("sz");

        // evaluate the nodal normals
        int N = m_surf.Nodes();
        vec3d nu(0, 0, 0);

        // loop over all elements to get average surface normal
        // (assumes that surface elements are all on same plane)
        for (int i = 0; i < m_surf.Elements(); ++i)
        {
            FESurfaceElement& el = m_surf.Element(i);
            nu += m_surf.SurfaceNormal(el, 0);
        }
        nu.unit();

        // create linear constraints
        // for a symmetry plane the constraint on (ux, uy, uz) is
        // nx*ux + ny*uy + nz*uz = 0
        for (int i = 0; i < N; ++i) {
            FENode node = m_surf.Node(i);
            if ((node.HasFlags(FENode::EXCLUDE) == false) && (node.m_rid == -1)) {
                FEAugLagLinearConstraint* pLC = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                for (int j = 0; j < 3; ++j) {
                    switch (j) {
                    case 0: pLC->AddDOF(node.GetID(), dofX, nu.x); break;
                    case 1: pLC->AddDOF(node.GetID(), dofY, nu.y); break;
                    case 2: pLC->AddDOF(node.GetID(), dofZ, nu.z); break;
                    default:
                        break;
                    }
                }
                // add the linear constraint to the system
                m_lc.add(pLC);
            }
        }

        // for nodes that belong to shells, also constraint the shell bottom face displacements
        for (int i = 0; i < N; ++i) {
            FENode node = m_surf.Node(i);
            if ((node.HasFlags(FENode::EXCLUDE) == false) && (node.HasFlags(FENode::SHELL)) && (node.m_rid == -1)) {
                FEAugLagLinearConstraint* pLC = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                for (int j = 0; j < 3; ++j) {
                    switch (j) {
                    case 0: pLC->AddDOF(node.GetID(), dofSX, nu.x); break;
                    case 1: pLC->AddDOF(node.GetID(), dofSY, nu.y); break;
                    case 2: pLC->AddDOF(node.GetID(), dofSZ, nu.z); break;
                    default:
                        break;
                    }
                }
                // add the linear constraint to the system
                m_lc.add(pLC);
            }
        }

        m_lc.Init();
        m_lc.Activate();

        m_binit = true;
    }
}

//-----------------------------------------------------------------------------
bool FESymmetryPlane::Init()
{
	// initialize surface
    return m_surf.Init();
}

//-----------------------------------------------------------------------------
void FESymmetryPlane::Serialize(DumpStream& ar) { m_lc.Serialize(ar); }
void FESymmetryPlane::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FESymmetryPlane::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FESymmetryPlane::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FESymmetryPlane::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
