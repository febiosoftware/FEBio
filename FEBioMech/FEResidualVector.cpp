#include "stdafx.h"
#include "FEResidualVector.h"
#include "FECore/FERigidSystem.h"
#include "FECore/FERigidBody.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"
#include <FECore/FELinearConstraintManager.h>
using namespace std;

//-----------------------------------------------------------------------------
FEResidualVector::FEResidualVector(FEModel& fem, vector<double>& R, vector<double>& Fr) : FEGlobalVector(fem, R, Fr)
{
}

//-----------------------------------------------------------------------------
FEResidualVector::~FEResidualVector()
{
}

//-----------------------------------------------------------------------------
void FEResidualVector::Assemble(vector<int>& en, vector<int>& elm, vector<double>& fe)
{
    
    FEModel& fem = GetFEModel();
    int dofX = fem.GetDOFIndex("x");
    int dofY = fem.GetDOFIndex("y");
    int dofZ = fem.GetDOFIndex("z");
    int dofU = fem.GetDOFIndex("u");
    int dofV = fem.GetDOFIndex("v");
    int dofW = fem.GetDOFIndex("w");
    
    vector<double>& R = m_R;
    
    int i, j, I, n, l;
    
    vec3d a, d;
    
    
    // get nodal DOFS
    //#pragma omp critical
    {
        // assemble the element residual into the global residual
        int ndof = fe.size();
        for (i=0; i<ndof; ++i)
        {
            
            I = elm[i];
            
            if ( I >= 0){
#pragma omp atomic
                R[I] += fe[i];
            }
            // TODO: Find another way to store reaction forces
            
            else if (-I-2 >= 0){
#pragma omp atomic
                m_Fr[-I-2] -= fe[i];
            }
        }
        
        
        int ndn = ndof / en.size();
        // if there are linear constraints we need to apply them
        
        
		// process linear constraints
		FELinearConstraintManager& LCM = m_fem.GetLinearConstraintManager();
		if (LCM.LinearConstraints())
		{
			LCM.AssembleResidual(R, en, elm, fe);
		}
        
        // If there are rigid bodies we need to look for rigid dofs
		FERigidSystem& rigid = *m_fem.GetRigidSystem();
        if (rigid.Objects() > 0)
        {
            int *lm;
            //#pragma omp critical
            for (i=0; i<ndof; i+=ndn)
            {
                FENode& node = m_fem.GetMesh().Node(en[i/ndn]);
                if (node.m_rid >= 0)
                {
                    
                    {
                        vec3d F(fe[i], fe[i+1], fe[i+2]);
                        
                        // this is an interface dof
                        // get the rigid body this node is connected to
                        FERigidBody& RB = *rigid.Object(node.m_rid);
                        lm = RB.m_LM;
                        
                        // add to total torque of this body
                        a = node.m_rt - RB.m_rt;
                        vec3d m = a ^ F;
                        vec3d f = F;
                        
                        if (node.HasFlags(FENode::SHELL)) {
                            vec3d d = node.m_d0 + node.get_vec3d(dofX, dofY, dofZ) - node.get_vec3d(dofU, dofV, dofW);
                            vec3d b = a - d;
                            vec3d Fd(fe[i+3], fe[i+4], fe[i+5]);
                            f += Fd;
                            m += b ^ Fd;
                        }
                        
                        n = lm[3];
                        if (n >= 0)
                        {
#pragma omp atomic
                            R[n] += m.x;
                        }
#pragma omp atomic
                        RB.m_Mr.x -= m.x;
                        n = lm[4];
                        if (n >= 0)
                        {
#pragma omp atomic
                            R[n] += m.y;
                        }
                        
#pragma omp atomic
                        RB.m_Mr.y -= m.y;
                        n = lm[5];
                        if (n >= 0)
                        {
#pragma omp atomic
                            R[n] += m.z;
                        }
#pragma omp atomic
                        RB.m_Mr.z -= m.z;
                        /*
                         // if the rotational degrees of freedom are constrained for a rigid node
                         // then we need to add an additional component to the residual
                         if (node.m_ID[m_dofRU] == lm[3])
                         {
                         d = node.m_Dt;
                         n = lm[3]; if (n >= 0) R[n] += d.y*F.z-d.z*F.y; RB.m_Mr.x -= d.y*F.z-d.z*F.y;
                         n = lm[4]; if (n >= 0) R[n] += d.z*F.x-d.x*F.z; RB.m_Mr.y -= d.z*F.x-d.x*F.z;
                         n = lm[5]; if (n >= 0) R[n] += d.x*F.y-d.y*F.x; RB.m_Mr.z -= d.x*F.y-d.y*F.x;
                         }
                         */
                        // add to global force vector
                        n = lm[0];
                        if (n >= 0)
                        {
#pragma omp atomic
                            R[n] += f.x;
                        }
#pragma omp atomic
                        RB.m_Fr.x -= f.x;
                        n = lm[1];
                        if (n >= 0) 
                        {
#pragma omp atomic
                            R[n] += f.y;
                        }
#pragma omp atomic
                        RB.m_Fr.y -= f.y;
                        
                        n = lm[2];
                        if (n >= 0)
                        {
#pragma omp atomic
                            R[n] += f.z;
                        }
#pragma omp atomic
                        RB.m_Fr.z -= f.z;
                    }
                }
            }
        }
    }
}
