#include "stdafx.h"
#include "FEBiphasicContactSurface.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FEBiphasicContactSurface::FEBiphasicContactSurface(FEModel* pfem) : FEContactSurface(pfem)
{
	m_dofP = -1;
}

//-----------------------------------------------------------------------------
FEBiphasicContactSurface::~FEBiphasicContactSurface()
{
}

//-----------------------------------------------------------------------------
bool FEBiphasicContactSurface::Init()
{
	// I want to use the FEModel class for this, but don't know how
	DOFS& dofs = GetFEModel()->GetDOFS();
	m_dofP = dofs.GetDOF("p");
	return FEContactSurface::Init();
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
double FEBiphasicContactSurface::GetFluidLoadSupport()
{
    double W = GetContactForce().norm();
    double Wp = GetFluidForce().norm();
    if (W == 0) return 0;
    return Wp/W;
}

//-----------------------------------------------------------------------------
void FEBiphasicContactSurface::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.assign(N*4, -1);

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
		if (m_dofP >= 0) lm[3*N+i] = id[m_dofP];
	}
}
