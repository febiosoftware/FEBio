#include "stdafx.h"
#include "FEResidualVector.h"
#include "FECore/FERigidBody.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"
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
        
        
        
        if (m_fem.m_LinC.size() > 0)
        {
		    DOFS& fedofs = m_fem.GetDOFS();
			int MAX_NDOFS = fedofs.GetTotalDOFS();

            // loop over all degrees of freedom of this element
            for (i=0; i<ndof; ++i)
            {
                // see if this dof belongs to a linear constraint
                n = MAX_NDOFS*(en[i/ndn]) + i%ndn;
                l = m_fem.m_LCT[n];
                
                if (l >= 0)
                {
                    // if so, get the linear constraint
                    FELinearConstraint& lc = *m_fem.m_LCA[l];
                    assert(elm[i] == -1);
                    
                    // now loop over all "slave" nodes and
                    // add the contribution to the residual
                    int ns = lc.slave.size();
                    list<FELinearConstraint::SlaveDOF>::iterator is = lc.slave.begin();
                    for (j=0; j<ns; ++j, ++is)
                    {
                        I = is->neq;
                        if (I >= 0)
                        {
                            double A = is->val;
#pragma omp atomic
                            R[I] += A*fe[i];
                        }
                    }
                }
            }
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
                        
                        n = lm[3];
                        if (n >= 0)
                        {
#pragma omp atomic
                            R[n] += a.y*F.z-a.z*F.y;
                        }
#pragma omp atomic
                        RB.m_Mr.x -= a.y*F.z-a.z*F.y;
                        n = lm[4];
                        if (n >= 0)
                        {
#pragma omp atomic
                            R[n] += a.z*F.x-a.x*F.z;
                        }
                        
#pragma omp atomic
                        RB.m_Mr.y -= a.z*F.x-a.x*F.z;
                        n = lm[5];
                        if (n >= 0)
                        {
#pragma omp atomic
                            R[n] += a.x*F.y-a.y*F.x;
                        }
#pragma omp atomic
                        RB.m_Mr.z -= a.x*F.y-a.y*F.x;
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
                            R[n] += F.x; 
                        }
#pragma omp atomic
                        RB.m_Fr.x -= F.x;
                        n = lm[1];
                        if (n >= 0) 
                        {
#pragma omp atomic
                            R[n] += F.y;
                        }
#pragma omp atomic
                        RB.m_Fr.y -= F.y;
                        
                        n = lm[2];
                        if (n >= 0)
                        {
#pragma omp atomic
                            R[n] += F.z; 
                        }
#pragma omp atomic
                        RB.m_Fr.z -= F.z;
                    }
                }
            }
        }
    }
}
