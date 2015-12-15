#include "stdafx.h"
#include "FEContactSurface.h"
#include <assert.h>

//-----------------------------------------------------------------------------
FEContactSurface::FEContactSurface(FEMesh* pm) : FESurface(pm) { m_pSibling = 0; }

//-----------------------------------------------------------------------------
FEContactSurface::~FEContactSurface() { m_pSibling = 0; }

//-----------------------------------------------------------------------------
void FEContactSurface::SetSibling(FEContactSurface* ps) { m_pSibling = ps; }

//-----------------------------------------------------------------------------
void FEContactSurface::GetNodalContactGap(int nface, double* pg) {}

//-----------------------------------------------------------------------------
void FEContactSurface::GetNodalContactPressure(int nface, double* pg) {}

//-----------------------------------------------------------------------------
void FEContactSurface::GetNodalContactTraction(int nface, vec3d* pt) {}

//-----------------------------------------------------------------------------
vec3d FEContactSurface::GetContactForce() { return vec3d(0,0,0); }

//-----------------------------------------------------------------------------
double FEContactSurface::GetContactArea() { return 0; }

//-----------------------------------------------------------------------------
void FEContactSurface::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*3);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

		lm[3*i  ] = id[DOF_X];
		lm[3*i+1] = id[DOF_Y];
		lm[3*i+2] = id[DOF_Z];
	}
}
