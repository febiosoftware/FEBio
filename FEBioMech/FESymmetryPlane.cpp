//
//  FESymmetryPlane.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 1/11/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#include "FESymmetryPlane.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FESymmetryPlane::FESymmetryPlane(FEModel* pfem) : FELinearConstraintSet(pfem), m_surf(&pfem->GetMesh())
{
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FESymmetryPlane::Activate()
{
    // don't forget to call base class
    FELinearConstraintSet::Activate();
}

//-----------------------------------------------------------------------------
bool FESymmetryPlane::Init()
{
    FEModel& fem = *FELinearConstraintSet::GetFEModel();
    DOFS& dofs = fem.GetDOFS();
    
    // initialize surface
    m_surf.Init();
    
    // evaluate the nodal normals
    int N = m_surf.Nodes();
    vec3d nu(0,0,0);
    
    // loop over all elements to get average surface normal
    // (assumes that surface elements are all on same plane)
    for (int i=0; i<m_surf.Elements(); ++i)
    {
        FESurfaceElement& el = m_surf.Element(i);
        nu += m_surf.SurfaceNormal(el, 0);
    }
    nu.unit();
    
    // create linear constraints
    // for a symmetry plane the constraint on (ux, uy, uz) is
    // nx*ux + ny*uy + nz*uz = 0
    for (int i=0; i<N; ++i) {
        FENode node = m_surf.Node(i);
        if (node.HasFlags(FENode::EXCLUDE) == false) {
            FEAugLagLinearConstraint* pLC = new FEAugLagLinearConstraint;
            for (int j=0; j<3; ++j) {
                FEAugLagLinearConstraint::DOF dof;
                dof.node = node.GetID() - 1;    // zero-based
                switch (j) {
                    case 0:
                        dof.bc = dofs.GetDOF("x");
                        dof.val = nu.x;
                        break;
                    case 1:
                        dof.bc = dofs.GetDOF("y");
                        dof.val = nu.y;
                        break;
                    case 2:
                        dof.bc = dofs.GetDOF("z");
                        dof.val = nu.z;
                        break;
                    default:
                        break;
                }
                pLC->m_dof.push_back(dof);
            }
            // add the linear constraint to the system
            add(pLC);
        }
    }
    
    // for nodes that belong to shells, also constraint the shell bottom face displacements
    for (int i=0; i<N; ++i) {
        FENode node = m_surf.Node(i);
        if ((node.HasFlags(FENode::EXCLUDE) == false) && (node.HasFlags(FENode::SHELL))) {
            FEAugLagLinearConstraint* pLC = new FEAugLagLinearConstraint;
            for (int j=0; j<3; ++j) {
                FEAugLagLinearConstraint::DOF dof;
                dof.node = node.GetID() - 1;    // zero-based
                switch (j) {
                    case 0:
                        dof.bc = dofs.GetDOF("u");
                        dof.val = nu.x;
                        break;
                    case 1:
                        dof.bc = dofs.GetDOF("v");
                        dof.val = nu.y;
                        break;
                    case 2:
                        dof.bc = dofs.GetDOF("w");
                        dof.val = nu.z;
                        break;
                    default:
                        break;
                }
                pLC->m_dof.push_back(dof);
            }
            // add the linear constraint to the system
            add(pLC);
        }
    }

    return true;
}
