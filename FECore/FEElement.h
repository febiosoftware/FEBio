// FEElement.h: interface for the FEElement class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEELEMENT_H__2EE38101_58E2_4FEB_B214_BB71B6FB15FB__INCLUDED_)
#define AFX_FEELEMENT_H__2EE38101_58E2_4FEB_B214_BB71B6FB15FB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEElementLibrary.h"
#include "FEElementTraits.h"
#include "FEMaterialPoint.h"
#include "FE_enum.h"
#include "FEException.h"

//-----------------------------------------------------------------------------
//! The FEElementState class stores the element state data. The state is defined
//! by a material point class for each of the integration points.
class FEElementState
{
public:
	//! default constructor
	FEElementState() {}

	//! destructor
	~FEElementState() { Clear(); }

	//! copy constructor
	FEElementState(const FEElementState& s);

	//! assignment operator
	FEElementState& operator = (const FEElementState& s);

	//! clear state data
	void Clear() { for (size_t i=0; i<m_data.size(); ++i) delete m_data[i]; m_data.clear(); }

	//! create 
	void Create(int n) { m_data.assign(n, static_cast<FEMaterialPoint*>(0) ); }

	//! operator for easy access to element data
	FEMaterialPoint*& operator [] (int n) { return m_data[n]; }

private:
	vector<FEMaterialPoint*>	m_data;
};

//-----------------------------------------------------------------------------
//! Base class for all element classes

//! From this class the different element classes are derived.

class FEElement
{
public:
	enum {MAX_NODES     = 20};	// max nr of nodes
	enum {MAX_INTPOINTS = 27};	// max nr of integration points

public:
	// default constructor
	FEElement() : m_pT(0) 
	{ 
		static int n = 1;
		m_nID = n++;
		m_nrigid = -1; 
	}

	virtual ~FEElement() {}

	bool IsRigid() { return (m_nrigid >= 0); }

	// Set the traits of an element
	virtual void SetTraits(FEElementTraits* ptraits)
	{
		m_pT = ptraits;
		m_node.resize(Nodes());
		m_State.Create(GaussPoints());
	}

	int GaussPoints() const { return m_pT->nint; } 
	int Nodes() const { return m_pT->neln; } 

	double* H(int n) { return m_pT->H[n]; }		// shape function values

	int Type() const { return m_pT->m_ntype; } 

	//! Get the element's material ID
	int GetMatID() const { return m_mat; } 

	//! Set the element's material ID
	void SetMatID(int id) { m_mat = id; }

	//! Set the type of the element
	void SetType(int ntype)
	{
		FEElementLibrary::SetElementTraits(*this, ntype);
	}

	void SetMaterialPointData(FEMaterialPoint* pmp, int n) { m_State[n] = pmp; }

	//! evaluate scalar field at integration point
	double Evaluate(double* fn, int n);

	//! evaluate scale field at integration point
	double Evaluate(vector<double>& fn, int n);

	//! evaluate vector field at integration point
	vec3d Evaluate(vec3d* vn, int n);

protected:
	int		m_mat;		//!< material index

public:

	int		m_nrigid;		//!< rigid body number that this element is attached to
	int		m_nID;			//!< element ID

	vector<int>			m_node;		//!< connectivity
	FEElementState		m_State;	//!< element state data
	FEElementTraits*	m_pT;		//!< pointer to element traits
};

//-----------------------------------------------------------------------------
//!  This class defines a solid element

class FESolidElement : public FEElement
{
public:
	//! default constructor
	FESolidElement(){}

	//! copy constructor
	FESolidElement(const FESolidElement& el);

	//! assignment operator
	FESolidElement& operator = (const FESolidElement& el);

	double* GaussWeights() { return &((FESolidElementTraits*)(m_pT))->gw[0]; }			// weights of integration points

	double* Gr(int n) { return ((FESolidElementTraits*)(m_pT))->Gr[n]; }	// shape function derivative to r
	double* Gs(int n) { return ((FESolidElementTraits*)(m_pT))->Gs[n]; }	// shape function derivative to s
	double* Gt(int n) { return ((FESolidElementTraits*)(m_pT))->Gt[n]; }	// shape function derivative to t

	double* Grr(int n) { return ((FESolidElementTraits*)(m_pT))->Grr[n]; }	// shape function 2nd derivative to rr
	double* Gsr(int n) { return ((FESolidElementTraits*)(m_pT))->Gsr[n]; }	// shape function 2nd derivative to sr
	double* Gtr(int n) { return ((FESolidElementTraits*)(m_pT))->Gtr[n]; }	// shape function 2nd derivative to tr
	
	double* Grs(int n) { return ((FESolidElementTraits*)(m_pT))->Grs[n]; }	// shape function 2nd derivative to rs
	double* Gss(int n) { return ((FESolidElementTraits*)(m_pT))->Gss[n]; }	// shape function 2nd derivative to ss
	double* Gts(int n) { return ((FESolidElementTraits*)(m_pT))->Gts[n]; }	// shape function 2nd derivative to ts
	
	double* Grt(int n) { return ((FESolidElementTraits*)(m_pT))->Grt[n]; }	// shape function 2nd derivative to rt
	double* Gst(int n) { return ((FESolidElementTraits*)(m_pT))->Gst[n]; }	// shape function 2nd derivative to st
	double* Gtt(int n) { return ((FESolidElementTraits*)(m_pT))->Gtt[n]; }	// shape function 2nd derivative to tt

	//! intialize element data
	void Init(bool bflag)
	{
		int nint = GaussPoints();
		for (int i=0; i<nint; ++i) m_State[i]->Init(bflag);
	}

	//! this function projects data from the gauss-points to the nodal points
	void project_to_nodes(double* ai, double* ao)
	{
		((FESolidElementTraits*)m_pT)->project_to_nodes(ai, ao);
	}
};

//-----------------------------------------------------------------------------
//!  This class defines a surface element

class FESurfaceElement : public FEElement
{
public:
	FESurfaceElement() { m_nelem = -1; m_lid = -1; }

	FESurfaceElement(const FESurfaceElement& el)
	{
		// set the traits of the element
		if (el.m_pT) SetTraits(el.m_pT);

		// copy base class data
		m_mat = el.m_mat;
		m_nrigid = el.m_nrigid;
		m_node = el.m_node;
		m_nID = el.m_nID;
		m_lid = el.m_lid;

		// copy surface element data
		m_nelem = el.m_nelem;
		m_lnode = el.m_lnode;
	}

	FESurfaceElement& operator = (const FESurfaceElement& el)
	{
		// make sure the element type is the same
		if (m_pT == 0) SetTraits(el.m_pT);
		else assert(m_pT == el.m_pT);

		// copy base class data
		m_mat = el.m_mat;
		m_nrigid = el.m_nrigid;
		m_node = el.m_node;
		m_nID = el.m_nID;
		m_lid = el.m_lid;

		// copy surface element data
		m_nelem = el.m_nelem;
		m_lnode = el.m_lnode;

		return (*this); 
	}

	virtual void SetTraits(FEElementTraits* pt)
	{
		// we don't allocate state data for surface elements
		m_pT = pt;
		m_node.resize(Nodes());
		m_lnode.resize(Nodes());
	}

	double* GaussWeights() { return &((FESurfaceElementTraits*)(m_pT))->gw[0]; }			// weights of integration points
	double gr(int n) { return ((FESurfaceElementTraits*)(m_pT))->gr[n]; }	// integration point coordinate r
	double gs(int n) { return ((FESurfaceElementTraits*)(m_pT))->gs[n]; }	// integration point coordinate  s

	double* Gr(int n) { return ((FESurfaceElementTraits*)(m_pT))->Gr[n]; }	// shape function derivative to r
	double* Gs(int n) { return ((FESurfaceElementTraits*)(m_pT))->Gs[n]; }	// shape function derivative to s

	double eval(double* d, int n)
	{
		double* N = H(n);
		int ne = Nodes();
		double a = 0;
		for (int i=0; i<ne; ++i) a += N[i]*d[i];
		return a;
	}

	double eval(double* d, double r, double s)
	{
		int n = Nodes();
		double H[FEElement::MAX_NODES];
		shape_fnc(H, r, s);
		double a = 0;
		for (int i=0; i<n; ++i) a += H[i]*d[i];
		return a;
	}

	vec3d eval(vec3d* d, double r, double s)
	{
		int n = Nodes();
		double H[FEElement::MAX_NODES];
		shape_fnc(H, r, s);
		vec3d a(0,0,0);
		for (int i=0; i<n; ++i) a += d[i]*H[i];
		return a;
	}

	vec3d eval(vec3d* d, int n)
	{
		int ne = Nodes();
		double* N = H(n);
		vec3d a(0,0,0);
		for (int i=0; i<ne; ++i) a += d[i]*N[i];
		return a;
	}

	double eval_deriv1(double* d, int j)
	{
		double* Hr = Gr(j);
		int n = Nodes();
		double s = 0;
		for (int i=0; i<n; ++i) s +=  Hr[i]*d[i];
		return s;
	}

	double eval_deriv2(double* d, int j)
	{
		double* Hs = Gs(j);
		int n = Nodes();
		double s = 0;
		for (int i=0; i<n; ++i) s +=  Hs[i]*d[i];
		return s;
	}

	double eval_deriv1(double* d, double r, double s)
	{
		double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
		shape_deriv(Hr, Hs, r, s);
		int n = Nodes();
		double a = 0;
		for (int i=0; i<n; ++i) a +=  Hr[i]*d[i];
		return a;
	}

	double eval_deriv2(double* d, double r, double s)
	{
		double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
		shape_deriv(Hr, Hs, r, s);
		int n = Nodes();
		double a = 0;
		for (int i=0; i<n; ++i) a +=  Hs[i]*d[i];
		return a;
	}

	void shape_fnc(double* H, double r, double s)
	{
		((FESurfaceElementTraits*)m_pT)->shape(H, r, s);
	}

	void shape_deriv(double* Gr, double* Gs, double r, double s)
	{
		((FESurfaceElementTraits*)m_pT)->shape_deriv(Gr, Gs, r, s);
	}

	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s)
	{
		((FESurfaceElementTraits*)m_pT)->shape_deriv2(Grr, Grs, Gss, r, s);
	}

	//! this function projects data from the gauss-points to the nodal points
	void project_to_nodes(double* ai, double* ao)
	{
		((FESurfaceElementTraits*)m_pT)->project_to_nodes(ai, ao);
	}

	bool HasNode(int n)
	{ 
		int l = Nodes(); 
		for (int i=0; i<l; ++i) 
			if (m_node[i] == n) return true; 
		return false;
	}

public:
	int		m_lid;			//!< local ID
	int		m_nelem;		//!< index of solid or shell element this surface element is a face of
	vector<int>	m_lnode;	//!< local node numbering (compared to m_node which is a global numbering)
};

//-----------------------------------------------------------------------------
//!  This class defines the shell element. 

//! A shell element is similar to a surface
//! element except that it has a thickness. 

class FEShellElement : public FEElement
{
public:
	FEShellElement(){}

	//! copy constructor
	FEShellElement(const FEShellElement& el);

	//! assignment operator
	FEShellElement& operator = (const FEShellElement& el);

	virtual void SetTraits(FEElementTraits* ptraits)
	{
		FEElement::SetTraits(ptraits);
		m_h0.resize(Nodes());
	}

	double* GaussWeights() { return &((FEShellElementTraits*)(m_pT))->gw[0]; }	// weights of integration points

	double* Hr(int n) { return ((FEShellElementTraits*)(m_pT))->Hr[n]; }	// shape function derivative to r
	double* Hs(int n) { return ((FEShellElementTraits*)(m_pT))->Hs[n]; }	// shape function derivative to s

	void Init(bool bflag)
	{
		int nint = GaussPoints();
		for (int i=0; i<nint; ++i) m_State[i]->Init(bflag);
	}

	double gr(int n) { return ((FEShellElementTraits*)(m_pT))->gr[n]; }
	double gs(int n) { return ((FEShellElementTraits*)(m_pT))->gs[n]; }
	double gt(int n) { return ((FEShellElementTraits*)(m_pT))->gt[n]; }

public:
	vector<double>	m_h0;	//!< initial shell thicknesses
};

//-----------------------------------------------------------------------------

class FETrussElement : public FEElement
{
public:
	FETrussElement(){}

	FETrussElement(const FETrussElement& el);

	FETrussElement& operator = (const FETrussElement& el);

	//! intialize element data
	void Init(bool bflag)
	{
		m_State[0]->Init(bflag);
	}

public:
	double	m_a0;	// cross-sectional area
};

//-----------------------------------------------------------------------------
//! Discrete element class

class FEDiscreteElement : public FEElement
{
public:
	FEDiscreteElement(){}
	FEDiscreteElement(const FEDiscreteElement& e);
	FEDiscreteElement& operator = (const FEDiscreteElement& e);
};

#endif // !defined(AFX_FEELEMENT_H__2EE38101_58E2_4FEB_B214_BB71B6FB15FB__INCLUDED_)
