/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



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
//! serialization
void FEBiphasicContactSurface::Serialize(DumpStream& ar)
{
	FEContactSurface::Serialize(ar);
	if (ar.IsShallow() == false) ar & m_dofP;
}

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
