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
    DOFS& dofs = *DOFS::GetInstance();
    int MAX_NDOFS = dofs.GetNDOFS();
    int MAX_CDOFS = dofs.GetCDOFS();
    
	int N = el.Nodes();
	lm.resize(N*(4+MAX_CDOFS));

	// get the DOF indices
	const int dof_x = dofs.GetDOF("x");
	const int dof_y = dofs.GetDOF("y");
	const int dof_z = dofs.GetDOF("z");
	const int dof_p = dofs.GetDOF("p");

	// pack the equation numbers
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];

		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[dof_x];
		lm[3*i+1] = id[dof_y];
		lm[3*i+2] = id[dof_z];

		// now the pressure dofs
		lm[3*N+i] = id[dof_p];

		// concentration dofs
		for (int k=0; k<MAX_CDOFS; ++k)
			lm[(4 + k)*N + i] = id[DOF_C+k];
	}
}
