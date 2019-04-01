#include "stdafx.h"
#include "FEConstraintNormalFlow.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FEConstraintNormalFlow::FEConstraintNormalFlow(FEModel* pfem) : FELinearConstraintSet(pfem), m_surf(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEConstraintNormalFlow::Activate()
{
    // don't forget to call base class
    FENLConstraint::Activate();
}

//-----------------------------------------------------------------------------
bool FEConstraintNormalFlow::Init()
{
    FEModel& fem = *GetFEModel();
    DOFS& dofs = fem.GetDOFS();
    
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
    
    // create linear constraints
    // for a surface with zero tangential velocity the constraints
    // on (vx, vy, vz) are
    //  (1-nx^2)*vx - nx*ny*vy - nx*nz*vz = 0
    // -nx*ny*vx + (1-ny^2)*vy - ny*nz*vz = 0
    // -nx*nz*vx - ny*nz*vy + (1-nz^2)*vz = 0
    for (int i=0; i<N; ++i) {

        //  (1-nx^2)*vx - nx*ny*vy - nx*nz*vz = 0
        FEAugLagLinearConstraint* pLC0 = new FEAugLagLinearConstraint;
        for (int j=0; j<3; ++j) {
            FEAugLagLinearConstraint::DOF dof;
            FENode node = m_surf.Node(i);
            dof.node = node.GetID() - 1;    // zero-based
            switch (j) {
                case 0:
                    dof.bc = dofs.GetDOF("wx");
                    dof.val = 1 - m_nn[i].x*m_nn[i].x;
                    break;
                case 1:
                    dof.bc = dofs.GetDOF("wy");
                    dof.val = -m_nn[i].x*m_nn[i].y;
                    break;
                case 2:
                    dof.bc = dofs.GetDOF("wz");
                    dof.val = -m_nn[i].x*m_nn[i].z;
                    break;
                default:
                    break;
            }
            pLC0->m_dof.push_back(dof);
        }
        // add the linear constraint to the system
        add(pLC0);
        
        // -nx*ny*vx + (1-ny^2)*vy - ny*nz*vz = 0
        FEAugLagLinearConstraint* pLC1 = new FEAugLagLinearConstraint;
        for (int j=0; j<3; ++j) {
            FEAugLagLinearConstraint::DOF dof;
            FENode node = m_surf.Node(i);
            dof.node = node.GetID() - 1;    // zero-based
            switch (j) {
                case 0:
                    dof.bc = dofs.GetDOF("wx");
                    dof.val = -m_nn[i].x*m_nn[i].y;
                    break;
                case 1:
                    dof.bc = dofs.GetDOF("wy");
                    dof.val = 1 - m_nn[i].y*m_nn[i].y;
                    break;
                case 2:
                    dof.bc = dofs.GetDOF("wz");
                    dof.val = -m_nn[i].y*m_nn[i].z;
                    break;
                default:
                    break;
            }
            pLC1->m_dof.push_back(dof);
        }
        // add the linear constraint to the system
        add(pLC1);
        
        // -nx*nz*vx - ny*nz*vy + (1-nz^2)*vz = 0
        FEAugLagLinearConstraint* pLC2 = new FEAugLagLinearConstraint;
        for (int j=0; j<3; ++j) {
            FEAugLagLinearConstraint::DOF dof;
            FENode node = m_surf.Node(i);
            dof.node = node.GetID() - 1;    // zero-based
            switch (j) {
                case 0:
                    dof.bc = dofs.GetDOF("wx");
                    dof.val = -m_nn[i].x*m_nn[i].z;
                    break;
                case 1:
                    dof.bc = dofs.GetDOF("wy");
                    dof.val = -m_nn[i].y*m_nn[i].z;
                    break;
                case 2:
                    dof.bc = dofs.GetDOF("wz");
                    dof.val = 1 - m_nn[i].z*m_nn[i].z;
                    break;
                default:
                    break;
            }
            pLC2->m_dof.push_back(dof);
        }
        // add the linear constraint to the system
        add(pLC2);
    }
    
    return true;
}

