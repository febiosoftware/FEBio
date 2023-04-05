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
#include "FEShellElement.h"
#include "DumpStream.h"
using namespace std;

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
	m_d0 = el.m_d0;
	m_g0[0] = el.m_g0[0]; m_g0[1] = el.m_g0[1]; m_g0[2] = el.m_g0[2];
	m_gt[0] = el.m_gt[0]; m_gt[1] = el.m_gt[1]; m_gt[2] = el.m_gt[2];
	m_gp[0] = el.m_gp[0]; m_gp[1] = el.m_gp[1]; m_gp[2] = el.m_gp[2];

	m_G0[0] = el.m_G0[0]; m_G0[1] = el.m_G0[1]; m_G0[2] = el.m_G0[2];
	m_Gt[0] = el.m_Gt[0]; m_Gt[1] = el.m_Gt[1]; m_Gt[2] = el.m_Gt[2];

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
	m_d0 = el.m_d0;
	m_g0[0] = el.m_g0[0]; m_g0[1] = el.m_g0[1]; m_g0[2] = el.m_g0[2];
	m_gt[0] = el.m_gt[0]; m_gt[1] = el.m_gt[1]; m_gt[2] = el.m_gt[2];
	m_gp[0] = el.m_gp[0]; m_gp[1] = el.m_gp[1]; m_gp[2] = el.m_gp[2];

	m_G0[0] = el.m_G0[0]; m_G0[1] = el.m_G0[1]; m_G0[2] = el.m_G0[2];
	m_Gt[0] = el.m_Gt[0]; m_Gt[1] = el.m_Gt[1]; m_Gt[2] = el.m_Gt[2];

	m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];

	return (*this);
}

void FEShellElement::SetTraits(FEElementTraits* ptraits)
{
	FEElement::SetTraits(ptraits);
	m_h0.assign(Nodes(), 0.0);
	m_ht.assign(Nodes(), 0.0);
	m_d0.assign(Nodes(), vec3d(0, 0, 0));
	m_g0[0].assign(GaussPoints(), vec3d(0, 0, 0));
	m_g0[1].assign(GaussPoints(), vec3d(0, 0, 0));
	m_g0[2].assign(GaussPoints(), vec3d(0, 0, 0));
	m_gt[0].assign(GaussPoints(), vec3d(0, 0, 0));
	m_gt[1].assign(GaussPoints(), vec3d(0, 0, 0));
	m_gt[2].assign(GaussPoints(), vec3d(0, 0, 0));
	m_gp[0].assign(GaussPoints(), vec3d(0, 0, 0));
	m_gp[1].assign(GaussPoints(), vec3d(0, 0, 0));
	m_gp[2].assign(GaussPoints(), vec3d(0, 0, 0));

	m_G0[0].assign(GaussPoints(), vec3d(0, 0, 0));
	m_G0[1].assign(GaussPoints(), vec3d(0, 0, 0));
	m_G0[2].assign(GaussPoints(), vec3d(0, 0, 0));
	m_Gt[0].assign(GaussPoints(), vec3d(0, 0, 0));
	m_Gt[1].assign(GaussPoints(), vec3d(0, 0, 0));
	m_Gt[2].assign(GaussPoints(), vec3d(0, 0, 0));
}

void FEShellElement::Serialize(DumpStream &ar)
{
	FEElement::Serialize(ar);
	if (ar.IsShallow() == false) {
		ar & m_h0;
		ar & m_d0;
		ar & m_g0[0] & m_g0[1] & m_g0[2];
		ar & m_gt[0] & m_gt[1] & m_gt[2];
		ar & m_gp[0] & m_gp[1] & m_gp[2];
		ar & m_G0[0] & m_G0[1] & m_G0[2];
		ar & m_Gt[0] & m_Gt[1] & m_Gt[2];
		ar & m_elem[0] & m_elem[1];
	}
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
	ar & m_fa;
	ar & m_Kaai;
	ar & m_alpha;
	ar & m_alphai;
	ar & m_alphat;
	ar & m_Kua;
	ar & m_Kwa;
	ar & m_E;
}
