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
FEElement::FEElement() : m_pT(0) 
{ 
	static int n = 1;
	m_nID = n++;
}

//-----------------------------------------------------------------------------
void FEElement::SetTraits(FEElementTraits* ptraits)
{
	m_pT = ptraits;
	m_node.resize(Nodes());
	m_lnode.resize(Nodes());
	m_State.Create(GaussPoints());
}

//-----------------------------------------------------------------------------
FESolidElement::FESolidElement(const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_node = el.m_node;
	m_lnode = el.m_lnode;
}

FESolidElement& FESolidElement::operator = (const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_node = el.m_node;
	m_lnode = el.m_lnode;

	return (*this);
}

//-----------------------------------------------------------------------------
FEShellElement::FEShellElement(const FEShellElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_node = el.m_node;
	m_lnode = el.m_lnode;

	// copy shell data
	m_h0 = el.m_h0;
	m_D0 = el.m_D0;
}

//! assignment operator
FEShellElement& FEShellElement::operator = (const FEShellElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_node = el.m_node;
	m_lnode = el.m_lnode;

	// copy shell data
	m_h0 = el.m_h0;
	m_D0 = el.m_D0;

	return (*this);
}

void FEShellElement::SetTraits(FEElementTraits* ptraits)
{
	FEElement::SetTraits(ptraits);
	m_h0.assign(Nodes(), 0.0);
	m_D0.resize(Nodes());
}

//-----------------------------------------------------------------------------
FEFergusonShellElement::FEFergusonShellElement(const FEFergusonShellElement& el)
{
    // set the traits of the element
    if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }
    
    // copy base class data
    m_mat = el.m_mat;
    m_nID = el.m_nID;
    m_node = el.m_node;
    m_lnode = el.m_lnode;
    
    // copy shell data
    m_h0 = el.m_h0;
    m_D0 = el.m_D0;
}

//! assignment operator
FEFergusonShellElement& FEFergusonShellElement::operator = (const FEFergusonShellElement& el)
{
    // set the traits of the element
    if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }
    
    // copy base class data
    m_mat = el.m_mat;
    m_nID = el.m_nID;
    m_node = el.m_node;
    m_lnode = el.m_lnode;
    
    // copy shell data
    m_h0 = el.m_h0;
    m_D0 = el.m_D0;
    
    return (*this);
}

void FEFergusonShellElement::SetTraits(FEElementTraits* ptraits)
{
    FEElement::SetTraits(ptraits);
    m_h0.assign(Nodes(), 0.0);
    m_D0.resize(Nodes());
}

//-----------------------------------------------------------------------------
FESurfaceElement::FESurfaceElement() 
{ 
	m_nelem = -1; m_lid = -1; 
}

FESurfaceElement::FESurfaceElement(const FESurfaceElement& el)
{
	// set the traits of the element
	if (el.m_pT) SetTraits(el.m_pT);

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_node = el.m_node;
	m_lnode = el.m_lnode;

	// copy surface element data
	m_lid = el.m_lid;
	m_nelem = el.m_nelem;
}

FESurfaceElement& FESurfaceElement::operator = (const FESurfaceElement& el)
{
	// make sure the element type is the same
	if (m_pT == 0) SetTraits(el.m_pT);
	else assert(m_pT == el.m_pT);

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_node = el.m_node;
	m_lnode = el.m_lnode;

	// copy surface element data
	m_lid = el.m_lid;
	m_nelem = el.m_nelem;

	return (*this); 
}

void FESurfaceElement::SetTraits(FEElementTraits* pt)
{
	// we don't allocate state data for surface elements
	m_pT = pt;
	m_node.resize(Nodes());
	m_lnode.resize(Nodes());
}

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
	m_node = el.m_node;
	m_lnode = el.m_lnode;

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
	m_node = el.m_node;
	m_lnode = el.m_lnode;

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
	m_node = el.m_node;
	m_lnode = el.m_lnode;
}

FEDiscreteElement& FEDiscreteElement::operator =(const FEDiscreteElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_node = el.m_node;
	m_lnode = el.m_lnode;

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
	m_node = el.m_node;
	m_lnode = el.m_lnode;
}

FEElement2D& FEElement2D::operator = (const FEElement2D& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nID = el.m_nID;
	m_node = el.m_node;
	m_lnode = el.m_lnode;

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
	m_node = el.m_node;
	m_lnode = el.m_lnode;
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
	m_node = el.m_node;
	m_lnode = el.m_lnode;

	return (*this);
}

void FELineElement::SetTraits(FEElementTraits* pt)
{
	// we don't allocate state data for surface elements
	m_pT = pt;
	m_node.resize(Nodes());
	m_lnode.resize(Nodes());
}
