#include "stdafx.h"
#include "FEFluidResidualVector.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"
#include <FECore/FELinearConstraintManager.h>
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
    
    int i, I;
    
    vec3d a, d;
    
    //#pragma omp critical
    {
        // assemble the element residual into the global residual
        int ndof = (int)fe.size();
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
        
        
        int ndn = ndof / (int)en.size();
        // if there are linear constraints we need to apply them
        

		// process linear constraints
		FELinearConstraintManager& LCM = m_fem.GetLinearConstraintManager();
		if (LCM.LinearConstraints())
		{
			LCM.AssembleResidual(R, en, elm, fe);
        }
    }
}
