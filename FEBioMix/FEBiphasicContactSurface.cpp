#include "stdafx.h"
#include "FEBiphasicContactSurface.h"
#include "FECore/DOFS.h"

//-----------------------------------------------------------------------------
FEBiphasicContactSurface::FEBiphasicContactSurface(FEMesh* pm) : FEContactSurface(pm)
{
}

//-----------------------------------------------------------------------------
FEBiphasicContactSurface::~FEBiphasicContactSurface()
{
}

//-----------------------------------------------------------------------------
void FEBiphasicContactSurface::GetNodalPressureGap(int nface, double* pg) { assert(false); }

//-----------------------------------------------------------------------------
vec3d FEBiphasicContactSurface::GetFluidForce()
{
	assert(false);
    return vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
void FEBiphasicContactSurface::UnpackLM(FEElement& el, vector<int>& lm)
{
    // get nodal DOFS
    DOFS& fedofs = *DOFS::GetInstance();
    int MAX_NDOFS = fedofs.GetNDOFS();
    int MAX_CDOFS = fedofs.GetCDOFS();
    
	int N = el.Nodes();
	lm.resize(N*(4+MAX_CDOFS));

	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];

		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[DOF_X];
		lm[3*i+1] = id[DOF_Y];
		lm[3*i+2] = id[DOF_Z];

		// now the pressure dofs
		lm[3*N+i] = id[DOF_P];

		// concentration dofs
		for (int k=0; k<MAX_CDOFS; ++k)
			lm[(4 + k)*N + i] = id[DOF_C+k];
	}
}
