#include "stdafx.h"
#include "FEMortarContactSurface.h"

//-----------------------------------------------------------------------------
FEMortarContactSurface::FEMortarContactSurface(FEMesh* pm) : FEContactSurface(pm)
{
}

//-----------------------------------------------------------------------------
bool FEMortarContactSurface::Init()
{
	if (FEContactSurface::Init() == false) return false;

	int NN = Nodes();
	m_gap.resize(NN, vec3d(0,0,0));
	return true;
}

//-----------------------------------------------------------------------------
void FEMortarContactSurface::UpdateNodalAreas()
{
	int NN = Nodes();
	int NF = Elements();
	m_A.resize(NN, 0.0);

	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = Element(i);
		double a = FaceArea(el);

		int nn = el.Nodes();
		double fa = a / (double) nn;
		for (int j=0; j<nn; ++j) m_A[el.m_lnode[j]] += fa;
	}

	for (int i=0; i<NN; ++i) m_A[i] = 1.0/m_A[i];
}
