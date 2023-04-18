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
#include "FESurfaceElement.h"
#include "DumpStream.h"

//=================================================================================================
// FESurfaceElement
//=================================================================================================

//-----------------------------------------------------------------------------
FESurfaceElement::FESurfaceElement()
{
	m_lid = -1;
	m_elem[0] = m_elem[1] = nullptr;
}

FESurfaceElement::FESurfaceElement(const FESurfaceElement& el) : FEElement(el)
{
	// set the traits of the element
	if (el.m_pT) SetTraits(el.m_pT);

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_lid = el.m_lid;
	m_node = el.m_node;
	m_lnode = el.m_lnode;
	m_lm = el.m_lm;
	m_val = el.m_val;
	m_status = el.m_status;

	// copy surface element data
	m_lid = el.m_lid;
	m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];
}

FESurfaceElement& FESurfaceElement::operator = (const FESurfaceElement& el)
{
	// make sure the element type is the same
	if (m_pT == 0) SetTraits(el.m_pT);
	else assert(m_pT == el.m_pT);

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_lid = el.m_lid;
	m_node = el.m_node;
	m_lnode = el.m_lnode;
	m_lm = el.m_lm;
	m_val = el.m_val;

	// copy surface element data
	m_lid = el.m_lid;
	m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];

	return (*this);
}

void FESurfaceElement::SetTraits(FEElementTraits* pt)
{
	m_pT = pt;
	m_node.resize(Nodes());
	m_lnode.resize(Nodes());
	m_State.Create(GaussPoints());
}

int FESurfaceElement::facet_edges() const
{
	int nn = Nodes(), nf = 0;
	switch (nn)
	{
	case 3:
	case 6:
	case 7:
		nf = 3;
		break;
	case 4:
	case 8:
	case 9:
		nf = 4;
		break;
	default:
		assert(false);
	}
	return nf;
}

void FESurfaceElement::facet_edge(int j, int* en) const
{
	int nn = Nodes();
	switch (nn)
	{
	case 3:
		en[0] = m_lnode[j];
		en[1] = m_lnode[(j + 1) % 3];
		break;
	case 6:
	case 7:
		en[0] = m_lnode[j];
		en[1] = m_lnode[j + 3];
		en[2] = m_lnode[(j + 1) % 3];
		break;
	case 4:
		en[0] = m_lnode[j];
		en[1] = m_lnode[(j + 1) % 4];
		break;
	case 8:
	case 9:
		en[0] = m_lnode[j];
		en[1] = m_lnode[j + 4];
		en[2] = m_lnode[(j + 1) % 4];
		break;
	}
}

double* FESurfaceElement::Gr(int order, int n) const { return (order >= 0 ? ((FESurfaceElementTraits*)(m_pT))->Gr_p[order][n] : ((FESurfaceElementTraits*)(m_pT))->Gr[n]); }	// shape function derivative to r
double* FESurfaceElement::Gs(int order, int n) const { return (order >= 0 ? ((FESurfaceElementTraits*)(m_pT))->Gs_p[order][n] : ((FESurfaceElementTraits*)(m_pT))->Gs[n]); }	// shape function derivative to s

void FESurfaceElement::Serialize(DumpStream& ar)
{
	FEElement::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_lid;
	// TODO: Serialize m_elem
}
