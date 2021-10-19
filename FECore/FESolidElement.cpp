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
#include "FESolidElement.h"
#include "DumpStream.h"

//=================================================================================================
// FESolidElement
//=================================================================================================

//-----------------------------------------------------------------------------
FESolidElement::FESolidElement(const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_lid = el.m_lid;
	m_node = el.m_node;
	m_lnode = el.m_lnode;
	m_lm = el.m_lm;
	m_val = el.m_val;
	m_status = el.m_status;
}

FESolidElement& FESolidElement::operator = (const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_lid = el.m_lid;
	m_node = el.m_node;
	m_lnode = el.m_lnode;
	m_lm = el.m_lm;
	m_val = el.m_val;
	m_status = el.m_status;

	return (*this);
}

void FESolidElement::SetTraits(FEElementTraits* pt)
{
	FEElement::SetTraits(pt);

	int ni = GaussPoints();
	m_J0i.resize(ni);
}

vec3d FESolidElement::evaluate(vec3d* v, double r, double s, double t) const
{
	double H[FEElement::MAX_NODES];
	shape_fnc(H, r, s, t);
	int neln = Nodes();
	vec3d p(0, 0, 0);
	for (int i = 0; i<neln; ++i) p += v[i] * H[i];
	return p;
}

double FESolidElement::evaluate(double* v, double r, double s, double t) const
{
	double H[FEElement::MAX_NODES];
	shape_fnc(H, r, s, t);
	int neln = Nodes();
	double p = 0.0;
	for (int i = 0; i<neln; ++i) p += v[i] * H[i];
	return p;
}

double* FESolidElement::Gr(int order, int n) const { return (order >= 0 ? ((FESolidElementTraits*)(m_pT))->m_Gr_p[order][n] : ((FESolidElementTraits*)(m_pT))->m_Gr[n]); }	// shape function derivative to r
double* FESolidElement::Gs(int order, int n) const { return (order >= 0 ? ((FESolidElementTraits*)(m_pT))->m_Gs_p[order][n] : ((FESolidElementTraits*)(m_pT))->m_Gs[n]); }	// shape function derivative to s
double* FESolidElement::Gt(int order, int n) const { return (order >= 0 ? ((FESolidElementTraits*)(m_pT))->m_Gt_p[order][n] : ((FESolidElementTraits*)(m_pT))->m_Gt[n]); }	// shape function derivative to t

void FESolidElement::Serialize(DumpStream& ar)
{
	FEElement::Serialize(ar);
	if (ar.IsShallow() == false)
	{
		ar & m_J0i;
		ar & m_bitfc;
	}
}
