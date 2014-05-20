// FEElement.cpp: implementation of the FEElement class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEElement.h"
#include <math.h>

//-----------------------------------------------------------------------------
FEElementState::FEElementState(const FEElementState& s)
{
	m_data.resize( s.m_data.size() );
	for (size_t i=0; i<m_data.size(); ++i) m_data[i] = s.m_data[i]->Copy();
}

FEElementState& FEElementState::operator = (const FEElementState& s)
{
	Clear();
	m_data.resize( s.m_data.size() );
	for (size_t i=0; i<m_data.size(); ++i) m_data[i] = s.m_data[i]->Copy();
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

vec3d FEElement::Evaluate(vec3d* vn, int n)
{
	double* Hn = H(n);
	vec3d v;
	const int N = Nodes();
	for (int i=0; i<N; ++i) v += vn[i]*Hn[i];
	return v;
}

//-----------------------------------------------------------------------------
FEElement::FEElement() : m_pT(0) 
{ 
	static int n = 1;
	m_nID = n++;
	m_nrigid = -1; 
}

//-----------------------------------------------------------------------------
void FEElement::SetTraits(FEElementTraits* ptraits)
{
	m_pT = ptraits;
	m_node.resize(Nodes());
	m_State.Create(GaussPoints());
}

//-----------------------------------------------------------------------------
FESolidElement::FESolidElement(const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
}

FESolidElement& FESolidElement::operator = (const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;

	return (*this);
}

//-----------------------------------------------------------------------------
FEShellElement::FEShellElement(const FEShellElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;

	// copy shell data
	m_h0 = el.m_h0;
}

//! assignment operator
FEShellElement& FEShellElement::operator = (const FEShellElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;

	// copy shell data
	m_h0 = el.m_h0;

	return (*this);
}

//-----------------------------------------------------------------------------
FETrussElement::FETrussElement(const FETrussElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;

	// truss data
	m_a0 = el.m_a0;
}

FETrussElement& FETrussElement::operator = (const FETrussElement& el) 
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;

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
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
}

FEDiscreteElement& FEDiscreteElement::operator =(const FEDiscreteElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;

	return (*this);
}
