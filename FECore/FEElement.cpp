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
#include "FEElement.h"
#include "DumpStream.h"
#include <math.h>

//-----------------------------------------------------------------------------
FEElementState::FEElementState(const FEElementState& s)
{
	m_data.resize( s.m_data.size() );
	for (size_t i=0; i<m_data.size(); ++i) 
	{
		if (s.m_data[i]) m_data[i] = s.m_data[i]->Copy(); else m_data[i] = 0;
	}
}

FEElementState& FEElementState::operator = (const FEElementState& s)
{
	Clear();
	m_data.resize( s.m_data.size() );
	for (size_t i=0; i<m_data.size(); ++i) 
	{
		if (s.m_data[i]) m_data[i] = s.m_data[i]->Copy(); else m_data[i] = 0;
	}
	return (*this);
}

//-----------------------------------------------------------------------------
double FEElement::Evaluate(double* fn, int n)
{
	double* Hn = H(n);
	double f = 0;
	const int N = Nodes();
	for (int i=0; i<N; ++i) f += Hn[i]*fn[i];
	return f;
}

double FEElement::Evaluate(vector<double>& fn, int n)
{
	double* Hn = H(n);
	double f = 0;
	const int N = Nodes();
	for (int i=0; i<N; ++i) f += Hn[i]*fn[i];
	return f;
}

vec2d FEElement::Evaluate(vec2d* vn, int n)
{
	double* Hn = H(n);
	vec2d v(0,0);
	const int N = Nodes();
	for (int i=0; i<N; ++i) v += vn[i]*Hn[i];
	return v;
}

vec3d FEElement::Evaluate(vec3d* vn, int n)
{
	double* Hn = H(n);
	vec3d v;
	const int N = Nodes();
	for (int i=0; i<N; ++i) v += vn[i]*Hn[i];
	return v;
}

//-----------------------------------------------------------------------------
double FEElement::Evaluate(double* fn, int order, int n)
{
	double* Hn = H(order, n);
	double f = 0;
	const int N = ShapeFunctions(order);
	for (int i = 0; i<N; ++i) f += Hn[i] * fn[i];
	return f;
}

double* FEElement::H(int order, int n)
{
	return m_pT->m_Hp[order][n];
}

int FEElement::ShapeFunctions(int order)
{
	return m_pT->ShapeFunctions(order);
}

//-----------------------------------------------------------------------------
bool FEElement::HasNode(int n) const
{
	int l = Nodes();
	for (int i = 0; i<l; ++i)
		if (m_node[i] == n) return true;
	return false;
}

//-----------------------------------------------------------------------------
int FEElement::FindNode(int n) const
{
	int l = Nodes();
	for (int i = 0; i<l; ++i)
		if (m_node[i] == n) return i;
	return -1;
}

//-----------------------------------------------------------------------------
FEElement::FEElement() : m_pT(0) 
{ 
	static int n = 1;
	m_nID = n++;
	m_lm = -1;
	m_val = 0.0;
	m_lid = -1;
	m_part = nullptr;
	m_status = ACTIVE;
}

//! get the element ID
int FEElement::GetID() const { return m_nID; }

//! set the element ID
void FEElement::SetID(int n) { m_nID = n; }

//! Get the element's material ID
int FEElement::GetMatID() const { return m_mat; }

//! Set the element's material ID
void FEElement::SetMatID(int id) { m_mat = id; }

//-----------------------------------------------------------------------------
void FEElement::SetTraits(FEElementTraits* ptraits)
{
	m_pT = ptraits;
	m_node.resize(Nodes());
	m_lnode.resize(Nodes());
	m_State.Create(GaussPoints());
}

//! serialize
void FEElement::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		int type = Type();
		ar << type;
		ar << m_nID << m_lid << m_mat;
		ar << m_node;
		ar << m_lnode;
		ar << m_lm << m_val;
		ar << m_status;
	}
	else
	{
		int ntype;
		ar >> ntype; SetType(ntype);
		ar >> m_nID >> m_lid >> m_mat;
		ar >> m_node;
		ar >> m_lnode;		
		ar >> m_lm >> m_val;
		ar >> m_status;
	}
}

//! return the nodes of the face
int FEElement::GetFace(int nface, int* nf) const
{
	int nn = -1;
	const int* en = &(m_node[0]);
	switch (Shape())
	{
	case ET_HEX8:
		nn = 4;
		switch (nface)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; break;
		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; break;
		case 3: nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; nf[3] = en[7]; break;
		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; break;
		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; break;
		}
		break;
	case ET_PENTA6:
		switch (nface)
		{
		case 0: nn = 4; nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[4]; nf[3] = en[3]; break;
		case 1: nn = 4; nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[5]; nf[3] = en[4]; break;
		case 2: nn = 4; nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[5]; nf[3] = en[2]; break;
		case 3: nn = 3; nf[0] = en[0]; nf[1] = en[2]; nf[2] = en[1]; nf[3] = en[1]; break;
		case 4: nn = 3; nf[0] = en[3]; nf[1] = en[4]; nf[2] = en[5]; nf[3] = en[5]; break;
		}
		break;
	case ET_PENTA15:
		switch (nface)
		{
		case 0: nn = 8; nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[4]; nf[3] = en[3]; nf[4] = en[6]; nf[5] = en[13]; nf[6] = en[9]; nf[7] = en[12]; break;
		case 1: nn = 8; nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[5]; nf[3] = en[4]; nf[4] = en[7]; nf[5] = en[14]; nf[6] = en[10]; nf[7] = en[13]; break;
		case 2: nn = 8; nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[5]; nf[3] = en[2]; nf[4] = en[12]; nf[5] = en[11]; nf[6] = en[14]; nf[7] = en[8]; break;
		case 3: nn = 6; nf[0] = en[0]; nf[1] = en[2]; nf[2] = en[1]; nf[3] = en[8]; nf[4] = en[7]; nf[5] = en[6]; break;
		case 4: nn = 6; nf[0] = en[3]; nf[1] = en[4]; nf[2] = en[5]; nf[3] = en[9]; nf[4] = en[10]; nf[5] = en[11]; break;
		}
		break;
	case ET_PYRA5:
		switch (nface)
		{
		case 0: nn = 3; nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[4]; break;
		case 1: nn = 3; nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[4]; break;
		case 2: nn = 3; nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[4]; break;
		case 3: nn = 3; nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; break;
		case 4: nn = 4; nf[0] = en[3]; nf[1] = en[2]; nf[2] = en[1]; nf[3] = en[0]; break;
		}
		break;
	case ET_TET4:
		nn = 3;
		switch (nface)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = nf[3] = en[3]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = nf[3] = en[3]; break;
		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = nf[3] = en[3]; break;
		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = nf[3] = en[0]; break;
		}
		break;
	case ET_TET10:
		nn = 6;
		switch (nface)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[3]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[7]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[3]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[8]; break;
		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = en[3]; nf[3] = en[6]; nf[4] = en[7]; nf[5] = en[9]; break;
		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = en[0]; nf[3] = en[5]; nf[4] = en[4]; nf[5] = en[6]; break;
		}
		break;
	case ET_TET15:
		nn = 7;
		switch (nface)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[3]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[7]; nf[6] = en[11]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[3]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[8]; nf[6] = en[12]; break;
		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = en[3]; nf[3] = en[6]; nf[4] = en[7]; nf[5] = en[9]; nf[6] = en[13]; break;
		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = en[0]; nf[3] = en[5]; nf[4] = en[4]; nf[5] = en[6]; nf[6] = en[10]; break;
		}
		break;
	case ET_TET20:
		nn = 10;
		switch (nface)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[3]; nf[3] = en[4]; nf[4] = en[5]; nf[5] = en[12]; nf[6] = en[13]; nf[7] = en[10]; nf[8] = en[11]; nf[9] = en[16]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[3]; nf[3] = en[6]; nf[4] = en[7]; nf[5] = en[14]; nf[6] = en[15]; nf[7] = en[13]; nf[8] = en[14]; nf[9] = en[17]; break;
		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = en[3]; nf[3] = en[9]; nf[4] = en[8]; nf[5] = en[10]; nf[6] = en[11]; nf[7] = en[14]; nf[8] = en[15]; nf[9] = en[18]; break;
		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = en[0]; nf[3] = en[7]; nf[4] = en[6]; nf[5] = en[5]; nf[6] = en[4]; nf[7] = en[10]; nf[8] = en[8]; nf[9] = en[19]; break;
		}
		break;
	case ET_HEX20:
		nn = 8;
		switch (nface)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[17]; nf[6] = en[12]; nf[7] = en[16]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[18]; nf[6] = en[13]; nf[7] = en[17]; break;
		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; nf[4] = en[10]; nf[5] = en[19]; nf[6] = en[14]; nf[7] = en[18]; break;
		case 3: nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; nf[3] = en[7]; nf[4] = en[11]; nf[5] = en[16]; nf[6] = en[15]; nf[7] = en[19]; break;
		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; nf[4] = en[11]; nf[5] = en[10]; nf[6] = en[9]; nf[7] = en[8]; break;
		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; nf[4] = en[12]; nf[5] = en[13]; nf[6] = en[14]; nf[7] = en[15]; break;
		}
		break;
	case ET_HEX27:
		nn = 9;
		switch (nface)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[17]; nf[6] = en[12]; nf[7] = en[16]; nf[8] = en[20]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[18]; nf[6] = en[13]; nf[7] = en[17]; nf[8] = en[21]; break;
		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; nf[4] = en[10]; nf[5] = en[19]; nf[6] = en[14]; nf[7] = en[18]; nf[8] = en[22]; break;
		case 3: nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; nf[3] = en[7]; nf[4] = en[11]; nf[5] = en[16]; nf[6] = en[15]; nf[7] = en[19]; nf[8] = en[23]; break;
		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; nf[4] = en[11]; nf[5] = en[10]; nf[6] = en[9]; nf[7] = en[8]; nf[8] = en[24]; break;
		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; nf[4] = en[12]; nf[5] = en[13]; nf[6] = en[14]; nf[7] = en[15]; nf[8] = en[25]; break;
		}
		break;
	case ET_QUAD4:
		nn = 4;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3];
		break;
	case ET_QUAD8:
		nn = 8;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5]; nf[6] = en[6]; nf[7] = en[7];
		break;
	case ET_QUAD9:
		nn = 9;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5]; nf[6] = en[6]; nf[7] = en[7]; nf[8] = en[8];
		break;
	case ET_TRI3:
		nn = 3;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2];
		break;
	case ET_TRI6:
		nn = 6;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5];
		break;
	case ET_TRI7:
		nn = 7;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5]; nf[6] = en[6];
		break;
	case ET_TRI10:
		nn = 7;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5]; nf[6] = en[6]; nf[7] = en[7]; nf[8] = en[8]; nf[9] = en[9];
		break;
	}

	return nn;
}

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
	vec3d p(0,0,0);
	for (int i=0; i<neln; ++i) p += v[i]*H[i];
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

void FESolidElement::Serialize(DumpStream& ar)
{
	FEElement::Serialize(ar);
	ar & m_J0i;
}

//=================================================================================================
// FEShellElement
//=================================================================================================

FEShellElement::FEShellElement()
{
	m_elem[0] = m_elem[1] = -1;
}

FEShellElement::FEShellElement(const FEShellElement& el)
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

	// copy shell data
	m_h0 = el.m_h0;
	m_ht = el.m_ht;

    m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];
}

//! assignment operator
FEShellElement& FEShellElement::operator = (const FEShellElement& el)
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

	// copy shell data
	m_h0 = el.m_h0;
	m_ht = el.m_ht;

	m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];

	return (*this);
}

void FEShellElement::SetTraits(FEElementTraits* ptraits)
{
	FEElement::SetTraits(ptraits);
	m_h0.assign(Nodes(), 0.0);
	m_ht.assign(Nodes(), 0.0);
}

void FEShellElement::Serialize(DumpStream &ar)
{
	FEElement::Serialize(ar);
	if (ar.IsShallow() == false) ar & m_h0;
	ar & m_ht;
}

//=================================================================================================
// FEShellElementOld
//=================================================================================================
FEShellElementOld::FEShellElementOld()
{
}

FEShellElementOld::FEShellElementOld(const FEShellElementOld& el) : FEShellElement(el)
{
	m_D0 = el.m_D0;
}

//! assignment operator
FEShellElementOld& FEShellElementOld::operator = (const FEShellElementOld& el)
{
	// copy base class
	FEShellElement::operator=(el);

	// copy this class data
	m_D0 = el.m_D0;

	return (*this);
}

void FEShellElementOld::SetTraits(FEElementTraits* ptraits)
{
	FEShellElement::SetTraits(ptraits);
	m_D0.resize(Nodes());
}

void FEShellElementOld::Serialize(DumpStream& ar)
{
	FEShellElement::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_D0;
}

//=================================================================================================
// FEShellElementNew
//=================================================================================================

FEShellElementNew::FEShellElementNew()
{
	
}

FEShellElementNew::FEShellElementNew(const FEShellElementNew& el) : FEShellElement(el)
{
	// TODO: What about all the EAS parameters?
}

//! assignment operator
FEShellElementNew& FEShellElementNew::operator = (const FEShellElementNew& el)
{
	FEShellElement::operator=(el);

	// TODO: What about all the EAS parameters?

	return (*this);
}

void FEShellElementNew::SetTraits(FEElementTraits* ptraits)
{
	FEShellElement::SetTraits(ptraits);

	// TODO: What about all the EAS parameters?
}

void FEShellElementNew::Serialize(DumpStream &ar)
{
	FEShellElement::Serialize(ar);

	if (ar.IsShallow())
	{
		if (ar.IsSaving()) {
			for (int i = 0; i<m_Kaai.rows(); ++i)
				for (int j = 0; j<m_Kaai.columns(); ++j)
					ar << m_Kaai(i, j);
			for (int i = 0; i<m_fa.rows(); ++i) ar << m_fa(i, 0);
			for (int i = 0; i<m_alpha.rows(); ++i) ar << m_alpha(i, 0);
			for (int i = 0; i<m_alphat.rows(); ++i) ar << m_alphat(i, 0);
			for (int i = 0; i<m_alphai.rows(); ++i) ar << m_alphai(i, 0);
			for (int k = 0; k<m_Kua.size(); ++k)
				for (int i = 0; i<m_Kua[k].rows(); ++i)
					for (int j = 0; j<m_Kua[k].columns(); ++j)
						ar << m_Kua[k](i, j);
			for (int k = 0; k<m_Kwa.size(); ++k)
				for (int i = 0; i<m_Kwa[k].rows(); ++i)
					for (int j = 0; j<m_Kwa[k].columns(); ++j)
						ar << m_Kwa[k](i, j);
			for (int k = 0; k<m_E.size(); ++k) ar << m_E[k];
		}
		else {
			for (int i = 0; i<m_Kaai.rows(); ++i)
				for (int j = 0; j<m_Kaai.columns(); ++j)
					ar >> m_Kaai(i, j);
			for (int i = 0; i<m_fa.rows(); ++i) ar >> m_fa(i, 0);
			for (int i = 0; i<m_alpha.rows(); ++i) ar >> m_alpha(i, 0);
			for (int i = 0; i<m_alphat.rows(); ++i) ar >> m_alphat(i, 0);
			for (int i = 0; i<m_alphai.rows(); ++i) ar >> m_alphai(i, 0);
			for (int k = 0; k<m_Kua.size(); ++k)
				for (int i = 0; i<m_Kua[k].rows(); ++i)
					for (int j = 0; j<m_Kua[k].columns(); ++j)
						ar >> m_Kua[k](i, j);
			for (int k = 0; k<m_Kwa.size(); ++k)
				for (int i = 0; i<m_Kwa[k].rows(); ++i)
					for (int j = 0; j<m_Kwa[k].columns(); ++j)
						ar >> m_Kwa[k](i, j);
			for (int k = 0; k<m_E.size(); ++k) ar >> m_E[k];
		}
	}
	else {
		if (ar.IsSaving()) {
			int nEAS = m_fa.rows();
			int neln = Nodes();
			int nint = GaussPoints();
			ar << nEAS;
			for (int i = 0; i<nEAS; ++i)
				for (int j = 0; j<nEAS; ++j)
					ar << m_Kaai(i, j);
			for (int i = 0; i<nEAS; ++i) ar << m_fa(i, 0);
			for (int i = 0; i<nEAS; ++i) ar << m_alpha(i, 0);
			for (int i = 0; i<nEAS; ++i) ar << m_alphat(i, 0);
			for (int i = 0; i<nEAS; ++i) ar << m_alphai(i, 0);
			for (int k = 0; k<neln; ++k)
				for (int i = 0; i<3; ++i)
					for (int j = 0; j<nEAS; ++j)
						ar << m_Kua[k](i, j);
			for (int k = 0; k<neln; ++k)
				for (int i = 0; i<3; ++i)
					for (int j = 0; j<nEAS; ++j)
						ar << m_Kwa[k](i, j);
			ar << m_E;
		}
		else {
			int nEAS;
			ar >> nEAS;
			int neln = Nodes();
			int nint = GaussPoints();
			for (int i = 0; i<nEAS; ++i)
				for (int j = 0; j<nEAS; ++j)
					ar >> m_Kaai(i, j);
			for (int i = 0; i<nEAS; ++i) ar >> m_fa(i, 0);
			for (int i = 0; i<nEAS; ++i) ar >> m_alpha(i, 0);
			for (int i = 0; i<nEAS; ++i) ar >> m_alphat(i, 0);
			for (int i = 0; i<nEAS; ++i) ar >> m_alphai(i, 0);
			for (int k = 0; k<neln; ++k)
				for (int i = 0; i<3; ++i)
					for (int j = 0; j<nEAS; ++j)
						ar >> m_Kua[k](i, j);
			for (int k = 0; k<neln; ++k)
				for (int i = 0; i<3; ++i)
					for (int j = 0; j<nEAS; ++j)
						ar >> m_Kwa[k](i, j);
			ar >> m_E;
		}
	}
}

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
            en[1] = m_lnode[(j+1)%3];
            break;
        case 6:
        case 7:
            en[0] = m_lnode[j];
            en[1] = m_lnode[(j+1)%3];
            en[2] = m_lnode[j+3];
            break;
        case 4:
            en[0] = m_lnode[j];
            en[1] = m_lnode[(j+1)%4];
            break;
        case 8:
        case 9:
            en[0] = m_lnode[j];
            en[1] = m_lnode[(j+1)%4];
            en[2] = m_lnode[j+4];
            break;
    }
}

void FESurfaceElement::Serialize(DumpStream& ar)
{
	FEElement::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_lid;
	// TODO: Serialize m_elem
}

//=================================================================================================
// FETrussElement
//=================================================================================================

//-----------------------------------------------------------------------------
FETrussElement::FETrussElement()
{
	m_a0 = 0.0;
}

FETrussElement::FETrussElement(const FETrussElement& el)
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

	// truss data
	m_a0 = el.m_a0;
}

FETrussElement& FETrussElement::operator = (const FETrussElement& el) 
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

	// copy truss data
	m_a0 = el.m_a0;

	return (*this); 
}

//-----------------------------------------------------------------------------
FEDiscreteElement::FEDiscreteElement(const FEDiscreteElement& el)
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
}

FEDiscreteElement& FEDiscreteElement::operator =(const FEDiscreteElement& el)
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

	return (*this);
}

//-----------------------------------------------------------------------------
FEElement2D::FEElement2D(const FEElement2D& el)
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
}

FEElement2D& FEElement2D::operator = (const FEElement2D& el)
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

	return (*this);
}

//-----------------------------------------------------------------------------
FELineElement::FELineElement()
{
	m_lid = -1;
}

FELineElement::FELineElement(const FELineElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy data
	m_lid = el.m_lid;

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_lid = el.m_lid;
	m_node = el.m_node;
	m_lnode = el.m_lnode;
	m_lm = el.m_lm;
	m_val = el.m_val;
}

FELineElement& FELineElement::operator = (const FELineElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy data
	m_lid = el.m_lid;

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_lid = el.m_lid;
	m_node = el.m_node;
	m_lnode = el.m_lnode;
	m_lm = el.m_lm;
	m_val = el.m_val;

	return (*this);
}

void FELineElement::SetTraits(FEElementTraits* pt)
{
	// we don't allocate state data for surface elements
	m_pT = pt;
	m_node.resize(Nodes());
	m_lnode.resize(Nodes());
}
