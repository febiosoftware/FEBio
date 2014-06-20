#include "stdafx.h"
#include "FEGlobalVector.h"
#include "FEModel.h"
#include "vec3d.h"

//-----------------------------------------------------------------------------
FEGlobalVector::FEGlobalVector(FEModel& fem, vector<double>& R, vector<double>& Fr) : m_fem(fem), m_R(R), m_Fr(Fr)
{
}

//-----------------------------------------------------------------------------
FEGlobalVector::~FEGlobalVector()
{

}

//-----------------------------------------------------------------------------
void FEGlobalVector::Assemble(vector<int>& en, vector<int>& elm, vector<double>& fe)
{
	vector<double>& R = m_R;

	// assemble the element residual into the global residual
	int ndof = fe.size();
	for (int i=0; i<ndof; ++i)
	{
		int I = elm[i];
		if ( I >= 0) R[I] += fe[i];
// TODO: Find another way to store reaction forces
		else if (-I-2 >= 0) m_Fr[-I-2] -= fe[i];
	}
}

//-----------------------------------------------------------------------------
void FEGlobalVector::AssembleRigid(int lm[6], double fe[6])
{
	vector<double>& R = m_R;
    
	int n;
    
    // add to global force vector
    n = lm[0]; if (n >= 0) R[n] += fe[0];
    n = lm[1]; if (n >= 0) R[n] += fe[1];
    n = lm[2]; if (n >= 0) R[n] += fe[2];
    
    // add to total torque of this body
    n = lm[3]; if (n >= 0) R[n] += fe[3];
    n = lm[4]; if (n >= 0) R[n] += fe[4];
    n = lm[5]; if (n >= 0) R[n] += fe[5];
}
