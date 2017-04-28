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

vec3d FESolidElement::evaluate(vec3d* v, double r, double s, double t) const
{
	double H[FEElement::MAX_NODES];
	shape_fnc(H, r, s, t);
	int neln = Nodes();
	vec3d p(0,0,0);
	for (int i=0; i<neln; ++i) p += v[i]*H[i];
	return p;
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
FESurfaceElement::FESurfaceElement() 
{ 
	m_lid = -1; 
	m_elem[0] = m_elem[1] = -1;
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
	m_node = el.m_node;
	m_lnode = el.m_lnode;

	// copy surface element data
	m_lid = el.m_lid;
	m_elem[0] = el.m_elem[0];
	m_elem[1] = el.m_elem[1];

	return (*this); 
}

void FESurfaceElement::SetTraits(FEElementTraits* pt)
{
	// we don't allocate state data for surface elements
	m_pT = pt;
	m_node.resize(Nodes());
	m_lnode.resize(Nodes());
}

void FESurfaceElement::project_to_nodes(vec3d* vi, vec3d* vo)
{
    int NELN = ((FESurfaceElementTraits*)m_pT)->neln;
    int NINT = ((FESurfaceElementTraits*)m_pT)->nint;
    double ai[FEElement::MAX_INTPOINTS], ao[FEElement::MAX_NODES];
    // along x
    for (int i=0; i<NINT; ++i) ai[i] = vi[i].x;
    ((FESurfaceElementTraits*)m_pT)->project_to_nodes(ai, ao);
    for (int j=0; j<NELN; ++j) vo[j].x = ao[j];
    // along y
    for (int i=0; i<NINT; ++i) ai[i] = vi[i].y;
    ((FESurfaceElementTraits*)m_pT)->project_to_nodes(ai, ao);
    for (int j=0; j<NELN; ++j) vo[j].y = ao[j];
    // along z
    for (int i=0; i<NINT; ++i) ai[i] = vi[i].z;
    ((FESurfaceElementTraits*)m_pT)->project_to_nodes(ai, ao);
    for (int j=0; j<NELN; ++j) vo[j].z = ao[j];
}

int FESurfaceElement::facet_edges()
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

void FESurfaceElement::facet_edge(int j, int* en)
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
