/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "FEPrescribedNormalDisplacement.h"
#include <FECore/FESurface.h>
#include <FECore/FENode.h>
#include "FEBioMech.h"

BEGIN_FECORE_CLASS(FEPrescribedNormalDisplacement, FEPrescribedSurface)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_hint, "surface_hint");
	ADD_PARAMETER(m_brelative, "relative");
END_FECORE_CLASS()

FEPrescribedNormalDisplacement::FEPrescribedNormalDisplacement(FEModel* fem) : FEPrescribedSurface(fem)
{
	m_scale = 0.0;
	m_hint = 0;

	// set the dof list
	// TODO: Can this be done in Init, since there is no error checking
	if (fem)
	{
		FEDofList dofs(fem);
		dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
		SetDOFList(dofs);
	}
}

// activation
void FEPrescribedNormalDisplacement::Activate()
{
	const FESurface& surf = *GetSurface();

	int N = surf.Nodes();
	m_node.resize(N);
	for (int i=0; i<N; ++i)
	{
		m_node[i].nodeId = surf.NodeIndex(i);
		m_node[i].normal = vec3d(0,0,0);
	}

	if (m_hint == 0)
	{
		int NF = surf.Elements();
		for (int i=0; i<NF; ++i)
		{
			const FESurfaceElement& el = surf.Element(i);
			int nn = el.Nodes();
			if ((nn==3) || (nn == 4))
			{
				for (int n=0; n<nn; ++n)
				{
					int i0 = el.m_lnode[n];
					int ip = el.m_lnode[(n+1)%nn];
					int im = el.m_lnode[(n + nn - 1)%nn];

					vec3d r0 = surf.Node(i0).m_r0;
					vec3d rp = surf.Node(ip).m_r0;
					vec3d rm = surf.Node(im).m_r0;

					vec3d nu = (rp - r0) ^ (rm - r0);

					m_node[i0].normal += nu;
				}
			}
			else if (nn == 6)
			{
				vec3d normals[6];

				// corner nodes
				for (int n = 0; n<3; ++n)
				{
					int i0 = el.m_lnode[n];
					int ip = el.m_lnode[(n + 1) % 3];
					int im = el.m_lnode[(n + 2) % 3];

					vec3d r0 = surf.Node(i0).m_r0;
					vec3d rp = surf.Node(ip).m_r0;
					vec3d rm = surf.Node(im).m_r0;

					vec3d nu = (rp - r0) ^ (rm - r0);

					normals[n] = nu;
				}

				// edge nodes
				for (int n=0; n<3; ++n)
				{
					int n0 = n + 3;
					int n1 = n;
					int n2 = (n+1) % 3;

					normals[n0] = (normals[n1] + normals[n2]) * 0.5;
				}

				for (int n=0; n<6; ++n) m_node[el.m_lnode[n]].normal += normals[n];
			}
		}

		for (int i = 0; i<N; ++i)
		{
			m_node[i].normal.unit();
		}
	}
	else
	{
		int NN = surf.Nodes();
		for (int i=0; i<NN; ++i)
		{
			vec3d ri = -surf.Node(i).m_r0;
			ri.unit();
			m_node[i].normal = ri;
		}
	}

	FEPrescribedSurface::Activate();
}

// return the values for node nodelid
void FEPrescribedNormalDisplacement::GetNodalValues(int nodelid, std::vector<double>& val)
{
	vec3d v = m_node[nodelid].normal*m_scale;
	val[0] = v.x;
	val[1] = v.y;
	val[2] = v.z;
}

// copy data from another class
void FEPrescribedNormalDisplacement::CopyFrom(FEBoundaryCondition* pbc)
{
	FEPrescribedNormalDisplacement* pnd = dynamic_cast<FEPrescribedNormalDisplacement*>(pbc);
	assert(pnd);
	m_scale = pnd->m_scale;
	m_node = pnd->m_node;
	CopyParameterListState(pnd->GetParameterList());
}
