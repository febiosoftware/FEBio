// FEElement.cpp: implementation of the FEElement class.
//
//////////////////////////////////////////////////////////////////////

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
		ar << Type();
		ar << m_nID << m_lid << m_mat;
		ar << m_node;
		ar << m_lnode;
		ar << m_lm << m_val;
	}
	else
	{
		int ntype;
		ar >> ntype; SetType(ntype);
		ar >> m_nID >> m_lid >> m_mat;
		ar >> m_node;
		ar >> m_lnode;		
		ar >> m_lm >> m_val;
	}
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

FESurfaceElement::FESurfaceElement(const FESurfaceElement& el)
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
