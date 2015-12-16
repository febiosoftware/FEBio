#include "stdafx.h"
#include "FEBiphasicContactSurface.h"
#include "FECore/DOFS.h"

//-----------------------------------------------------------------------------
FEBiphasicContactSurface::FEBiphasicContactSurface(FEMesh* pm) : FEContactSurface(pm)
{
	// I want to use the FEModel class for this, but don't know how
	DOFS& dofs = *DOFS::GetInstance();
	m_dofP = dofs.GetDOF("p");
	m_dofC = dofs.GetDOF("c");
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
    int MAX_CDOFS = dofs.GetDOFSize("c");
    
	int N = el.Nodes();
	lm.resize(N*(4+MAX_CDOFS));

	// pack the equation numbers
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];

		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[m_dofX];
		lm[3*i+1] = id[m_dofY];
		lm[3*i+2] = id[m_dofZ];

		// now the pressure dofs
		lm[3*N+i] = id[m_dofP];

		// concentration dofs
		for (int k=0; k<MAX_CDOFS; ++k)
			lm[(4 + k)*N + i] = id[m_dofC+k];
	}
}
