#include "stdafx.h"
#include "FEContactSurface.h"
#include "FECore/FEModel.h"
#include <assert.h>

//-----------------------------------------------------------------------------
FEContactSurface::FEContactSurface(FEModel* pfem) : FESurface(&pfem->GetMesh()), m_pfem(pfem)
{
	m_pSibling = 0; 
	m_dofX = -1;
	m_dofY = -1;
	m_dofZ = -1;
}

//-----------------------------------------------------------------------------
FEContactSurface::~FEContactSurface() { m_pSibling = 0; }

//-----------------------------------------------------------------------------
bool FEContactSurface::Init()
{
	// I want to use the FEModel class for this, but don't know how
	DOFS& dofs = GetFEModel()->GetDOFS();
	m_dofX = dofs.GetDOF("x");
	m_dofY = dofs.GetDOF("y");
	m_dofZ = dofs.GetDOF("z");

	return FESurface::Init();
}

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

		lm[3*i  ] = id[m_dofX];
		lm[3*i+1] = id[m_dofY];
		lm[3*i+2] = id[m_dofZ];
	}
}
