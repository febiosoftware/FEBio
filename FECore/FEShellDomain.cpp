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
#include "FEShellDomain.h"
#include "FEMesh.h"
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! constructor
FEShellDomain::FEShellDomain(FEModel* fem) : FEDomain(FE_DOMAIN_SHELL, fem)
{
}

//-----------------------------------------------------------------------------
void FEShellDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	ForEachMaterialPoint([&](FEMaterialPoint& mp) {
		mp.Update(timeInfo);
	});
}

//-----------------------------------------------------------------------------
void FEShellDomain::Reset()
{
	ForEachShellElement([](FEShellElement& el) {
		int ni = el.GaussPoints();
		for (int j = 0; j<ni; ++j) el.GetMaterialPoint(j)->Init();

		int ne = el.Nodes();
		for (int j = 0; j<ne; ++j) el.m_ht[j] = el.m_h0[j];
	});
}

//-----------------------------------------------------------------------------
void FEShellDomain::InitShells()
{
	ForEachShellElement([](FEShellElement& el) {
		int n = el.Nodes();
		for (int j = 0; j<n; ++j) el.m_ht[j] = el.m_h0[j];
	});
}

//-----------------------------------------------------------------------------
void FEShellDomain::ForEachShellElement(std::function<void(FEShellElement& el)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i) f(Element(i));
}

//=================================================================================================

FEShellDomainOld::FEShellDomainOld(FEModel* fem) : FEShellDomain(fem)
{
}

//-----------------------------------------------------------------------------
void FEShellDomainOld::Create(int nelems, int elemType)
{
	m_Elem.resize(nelems);
	for (int i = 0; i < nelems; ++i)
	{
		FEShellElementOld& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	if (elemType != -1)
		for (int i=0; i<nelems; ++i) m_Elem[i].SetType(elemType);
}

//-----------------------------------------------------------------------------
double FEShellDomainOld::Volume(FEShellElement& se)
{
	FEShellElementOld& el = static_cast<FEShellElementOld&>(se);

	int neln = el.Nodes();

	// initial nodal coordinates and directors
	vec3d r0[FEElement::MAX_NODES], D0[FEElement::MAX_NODES];
	for (int i = 0; i<neln; ++i)
	{
		r0[i] = Node(el.m_lnode[i]).m_r0;
		D0[i] = el.m_D0[i];
	}

	int nint = el.GaussPoints();
	double *w = el.GaussWeights();
	double V = 0;
	vec3d g[3];
	for (int n = 0; n<nint; ++n)
	{
		// jacobian matrix
		double eta = el.gt(n);

		double* Mr = el.Hr(n);
		double* Ms = el.Hs(n);
		double* M = el.H(n);

		// evaluate covariant basis vectors
		g[0] = g[1] = g[2] = vec3d(0, 0, 0);
		for (int i = 0; i<neln; ++i)
		{
			g[0] += (r0[i] + D0[i] * eta / 2)*Mr[i];
			g[1] += (r0[i] + D0[i] * eta / 2)*Ms[i];
			g[2] += D0[i] * (M[i] / 2);
		}

		mat3d J = mat3d(g[0].x, g[1].x, g[2].x,
			g[0].y, g[1].y, g[2].y,
			g[0].z, g[1].z, g[2].z);

		// calculate the determinant
		double detJ0 = J.det();

		V += detJ0*w[n];
	}

	return V;
}

//-----------------------------------------------------------------------------
//! Calculate all shell normals (i.e. the shell directors).
//! And find shell nodes
void FEShellDomainOld::InitShells()
{
	FEShellDomain::InitShells();

	FEMesh& mesh = *GetMesh();
	for (int i = 0; i<Elements(); ++i)
	{
		FEShellElementOld& el = ShellElement(i);
		int ne = el.Nodes();
		for (int j = 0; j<ne; ++j)
		{
			vec3d d0 = mesh.Node(el.m_node[j]).m_d0;
			d0.unit();
			el.m_D0[j] = d0 * el.m_h0[j];
		}
	}
}

//=================================================================================================

FEShellDomainNew::FEShellDomainNew(FEModel* fem) : FEShellDomain(fem)
{
}

//-----------------------------------------------------------------------------
void FEShellDomainNew::Create(int nelems, int elemType)
{
	m_Elem.resize(nelems);
	for (int i = 0; i < nelems; ++i)
	{
		FEShellElementNew& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	if (elemType != -1)
		for (int i = 0; i<nelems; ++i) m_Elem[i].SetType(elemType);
}

//-----------------------------------------------------------------------------
double FEShellDomainNew::Volume(FEShellElement& se)
{
	FEShellElementNew& el = static_cast<FEShellElementNew&>(se);

	int neln = el.Nodes();

	// initial nodal coordinates and directors
	vec3d r0[FEElement::MAX_NODES], D0[FEElement::MAX_NODES];
	for (int i = 0; i<neln; ++i)
	{
		r0[i] = Node(el.m_lnode[i]).m_r0;
		D0[i] = Node(el.m_lnode[i]).m_d0;
	}

	int nint = el.GaussPoints();
	double *w = el.GaussWeights();
	double V = 0;
	vec3d g[3];
	for (int n = 0; n<nint; ++n)
	{
		// jacobian matrix
		double eta = el.gt(n);

		double* Mr = el.Hr(n);
		double* Ms = el.Hs(n);
		double* M = el.H(n);

		// evaluate covariant basis vectors
		g[0] = g[1] = g[2] = vec3d(0, 0, 0);
		for (int i = 0; i<neln; ++i)
		{
			g[0] += (r0[i] + D0[i] * eta / 2)*Mr[i];
			g[1] += (r0[i] + D0[i] * eta / 2)*Ms[i];
			g[2] += D0[i] * (M[i] / 2);
		}

		mat3d J = mat3d(g[0].x, g[1].x, g[2].x,
			g[0].y, g[1].y, g[2].y,
			g[0].z, g[1].z, g[2].z);

		// calculate the determinant
		double detJ0 = J.det();

		V += detJ0*w[n];
	}

	return V;
}


