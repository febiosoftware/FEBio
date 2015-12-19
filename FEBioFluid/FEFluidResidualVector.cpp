#include "stdafx.h"
#include "FEFluidResidualVector.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"
using namespace std;

//-----------------------------------------------------------------------------
FEFluidResidualVector::FEFluidResidualVector(FEModel& fem, vector<double>& R, vector<double>& Fr) : FEGlobalVector(fem, R, Fr)
{
}

//-----------------------------------------------------------------------------
FEFluidResidualVector::~FEFluidResidualVector()
{
}

//-----------------------------------------------------------------------------
void FEFluidResidualVector::Assemble(vector<int>& en, vector<int>& elm, vector<double>& fe)
{
    
    vector<double>& R = m_R;
    
    int i, j, I, n, l;
    
    vec3d a, d;
    
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
		    // get nodal DOFS
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
    }
}
