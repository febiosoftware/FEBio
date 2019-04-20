/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEFluidTractionLoad.h"
#include "FECore/FEModel.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidTractionLoad, FESurfaceLoad)
	ADD_PARAMETER(m_scale, "scale"   );
	ADD_PARAMETER(m_TC   , "traction");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidTractionLoad::FEFluidTractionLoad(FEModel* pfem) : FESurfaceLoad(pfem), m_TC(FE_VEC3D)
{
	m_scale = 1.0;

	m_dofWX = pfem->GetDOFIndex("wx");
	m_dofWY = pfem->GetDOFIndex("wy");
	m_dofWZ = pfem->GetDOFIndex("wz");
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidTractionLoad::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	m_TC.Create(ps); 
}

//-----------------------------------------------------------------------------
void FEFluidTractionLoad::UnpackLM(FEElement& el, vector<int>& lm)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	int N = el.Nodes();
	lm.resize(N*3);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = mesh.Node(n);
		vector<int>& id = node.m_ID;

		lm[3*i  ] = id[m_dofWX];
		lm[3*i+1] = id[m_dofWY];
		lm[3*i+2] = id[m_dofWZ];
	}
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
void FEFluidTractionLoad::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
	vector<double> fe;
	vector<int> elm;

	vec3d r0[FEElement::MAX_NODES];
	vec3d tn[FEElement::MAX_NODES];

	int i, n;
	int N = m_psurf->Elements();
	for (int iel=0; iel<N; ++iel)
	{
		FESurfaceElement& el = m_psurf->Element(iel);

        int ndof = 3*el.Nodes();
		fe.resize(ndof);

		// nr integration points
		int nint = el.GaussPoints();

		// nr of element nodes
		int neln = el.Nodes();

		// nodal coordinates
		for (i=0; i<neln; ++i)
		{
			r0[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_r0;
			tn[i] = m_TC.value<vec3d>(iel, i)*m_scale;
		}

		double* Gr, *Gs;
		double* N;
		double* w  = el.GaussWeights();

		vec3d dxr, dxs;

		// repeat over integration points
		zero(fe);
		for (n=0; n<nint; ++n)
		{
			N  = el.H(n);
			Gr = el.Gr(n);
			Gs = el.Gs(n);

			// calculate the tangent vectors
			dxr = dxs = vec3d(0,0,0);
			vec3d t(0,0,0);
			for (i=0; i<neln; ++i) 
			{
				dxr.x += Gr[i]*r0[i].x;
				dxr.y += Gr[i]*r0[i].y;
				dxr.z += Gr[i]*r0[i].z;

				dxs.x += Gs[i]*r0[i].x;
				dxs.y += Gs[i]*r0[i].y;
				dxs.z += Gs[i]*r0[i].z;

				t += tn[i]*N[i];
			}

			vec3d f = t*((dxr ^ dxs).norm()*w[n]);

			for (i=0; i<neln; ++i)
			{
                fe[3*i  ] += N[i]*f.x;
                fe[3*i+1] += N[i]*f.y;
                fe[3*i+2] += N[i]*f.z;
			}
		}

		// get the element's LM vector and adjust it
		UnpackLM(el, elm);
        
		// add element force vector to global force vector
		R.Assemble(el.m_node, elm, fe);
	}
}
