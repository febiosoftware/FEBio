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
//! get the current nodal coordinates
void FEShellDomain::GetCurrentNodalCoordinates(const FEShellElement& el, vec3d* rt, const bool back)
{
    int neln = el.Nodes();
    if (!back)
        for (int i = 0; i<neln; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
    else
        for (int i = 0; i<neln; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).st();
}

//-----------------------------------------------------------------------------
//! get the current nodal coordinates
void FEShellDomain::GetCurrentNodalCoordinates(const FEShellElement& el, vec3d* rt, double alpha, const bool back)
{
    int neln = el.Nodes();
    if (!back) {
        for (int i = 0; i<neln; ++i) {
            FENode& nd = m_pMesh->Node(el.m_node[i]);
            rt[i] = nd.m_rt*alpha + nd.m_rp*(1 - alpha);
        }
    }
    else {
        for (int i = 0; i<neln; ++i) {
            FENode& nd = m_pMesh->Node(el.m_node[i]);
            rt[i] = nd.st()*alpha + nd.sp()*(1 - alpha);
        }
    }
}

//-----------------------------------------------------------------------------
//! get the reference nodal coordinates
void FEShellDomain::GetReferenceNodalCoordinates(const FEShellElement& el, vec3d* r0, const bool back)
{
    int neln = el.Nodes();
    if (!back)
        for (int i = 0; i<neln; ++i) r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
    else
        for (int i = 0; i<neln; ++i) r0[i] = m_pMesh->Node(el.m_node[i]).s0();
}

//-----------------------------------------------------------------------------
//! get the previous nodal coordinates
void FEShellDomain::GetPreviousNodalCoordinates(const FEShellElement& el, vec3d* rp, const bool back)
{
    int neln = el.Nodes();
    if (!back)
        for (int i = 0; i<neln; ++i) rp[i] = m_pMesh->Node(el.m_node[i]).m_rp;
    else
        for (int i = 0; i<neln; ++i) rp[i] = m_pMesh->Node(el.m_node[i]).sp();
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
bool FEShellDomainOld::Create(int nelems, FE_Element_Spec espec)
{
	m_Elem.resize(nelems);
	for (int i = 0; i < nelems; ++i)
	{
		FEShellElementOld& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	if (espec.etype != FE_ELEM_INVALID_TYPE)
		for (int i=0; i<nelems; ++i) m_Elem[i].SetType(espec.etype);

	return true;
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

BEGIN_FECORE_CLASS(FEShellDomainNew, FEShellDomain)
	ADD_PARAMETER(m_h0, "shell_thickness");
END_FECORE_CLASS();

FEShellDomainNew::FEShellDomainNew(FEModel* fem) : FEShellDomain(fem)
{
	m_h0 = 0.0;
}

//-----------------------------------------------------------------------------
bool FEShellDomainNew::Create(int nelems, FE_Element_Spec espec)
{
	m_Elem.resize(nelems);
	for (int i = 0; i < nelems; ++i)
	{
		FEShellElementNew& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	if (espec.etype != FE_ELEM_INVALID_TYPE)
		for (int i = 0; i<nelems; ++i) m_Elem[i].SetType(espec.etype);

	return true;
}

//-----------------------------------------------------------------------------
void FEShellDomainNew::AssignDefaultShellThickness()
{
	double h0 = DefaultShellThickness();
	if (h0 <= 0.0) return;

	for (int j = 0; j < Elements(); ++j)
	{
		FEShellElement& el = Element(j);
		int ne = el.Nodes();
		for (int n = 0; n < ne; ++n) el.m_ht[n] = el.m_h0[n] = h0;
	}
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


