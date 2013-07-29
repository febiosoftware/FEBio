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
