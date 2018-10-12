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

class FEMesh;
//-----------------------------------------------------------------------------
class FEElementTraits;
class FEDomain;

//-----------------------------------------------------------------------------
//! The FEElementState class stores the element state data. The state is defined
//! by a material point class for each of the integration points.
class FECORE_API FEElementState
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

class FECORE_API FEElement
{
public:
	enum {MAX_NODES     = 27};	// max nr of nodes
	enum {MAX_INTPOINTS = 27};	// max nr of integration points

public:
	//! default constructor
	FEElement();

	//! destructor
	virtual ~FEElement() {}

	//! get the element ID
	int GetID() const { return m_nID; }

	//! set the element ID
	void SetID(int n) { m_nID = n; }

	//! Get the element's material ID
	int GetMatID() const { return m_mat; }

	//! Set the element's material ID
	void SetMatID(int id) { m_mat = id; }

	//Get the domain that contains this element
	FEDomain * GetDomain() const { return m_dom; }

	//Set the domain that contains this element
	void SetDomain(FEDomain * dom){ m_dom = dom; }

	//! Set the Local ID
	void SetLocalID(int lid) { m_lid = lid; }

	//! Get the local ID
	int GetLocalID() const { return m_lid; }

public:
	//! Set the type of the element
	void SetType(int ntype) { FEElementLibrary::SetElementTraits(*this, ntype); }

	//! Set the traits of an element
	virtual void SetTraits(FEElementTraits* ptraits);

	//! Get the element traits
	FEElementTraits* GetTraits() { return m_pT; }

	//! return number of nodes
	int Nodes() const { return m_pT->neln; } 

	//! return the element class
	int Class() const { return m_pT->Class(); }

	//! return the element shape
	int Shape() const { return m_pT->Shape(); }

	//! return the type of element
	int Type() const { return m_pT->Type(); }

	//! return number of integration points
	int GaussPoints() const { return m_pT->nint; } 

	//! shape function values
	double* H(int n) { return m_pT->H[n]; }

public:
	//! Get the material point data
	FEMaterialPoint* GetMaterialPoint(int n) { return m_State[n]; }

	//! set the material point data
	void SetMaterialPointData(FEMaterialPoint* pmp, int n)
	{ 
		pmp->m_elem = this;
		pmp->m_index = n;
		m_State[n] = pmp; 
	}

	//! serialize
	//! NOTE: state data is not serialized by the element. This has to be done by the domains.
	virtual void Serialize(DumpStream& ar);

public:
	//! evaluate scalar field at integration point
	double Evaluate(double* fn, int n);

	//! evaluate scale field at integration point
	double Evaluate(vector<double>& fn, int n);

	//! evaluate vector field at integration point
	vec2d Evaluate(vec2d* vn, int n);

	//! evaluate vector field at integration point
	vec3d Evaluate(vec3d* vn, int n);

	// see if this element has the node n
    bool HasNode(int n) const;

	// find local element index of node n    
    int FindNode(int n) const;

	// project data to nodes
	void project_to_nodes(double* ai, double* ao) const { m_pT->project_to_nodes(ai, ao); }
	void project_to_nodes(vec3d*  ai, vec3d*  ao) const { m_pT->project_to_nodes(ai, ao); }
	void project_to_nodes(mat3ds* ai, mat3ds* ao) const { m_pT->project_to_nodes(ai, ao); }
	void project_to_nodes(mat3d*  ai, mat3d*  ao) const { m_pT->project_to_nodes(ai, ao); }

protected:
	int		m_nID;		//!< element ID
	int		m_lid;		//!< local ID
	int		m_mat;		//!< material index
	FEDomain * m_dom;	//!< parent domain

public:
	vector<int>		m_node;		//!< connectivity
	int		m_lm;

	// This array stores the local node numbers, that is the node numbers
	// into the node list of a domain.
	vector<int>		m_lnode;	//!< local connectivity

protected:
	FEElementState		m_State;	//!< element state data
	FEElementTraits*	m_pT;		//!< pointer to element traits
};

//-----------------------------------------------------------------------------
//!  This class defines a solid element

class FECORE_API FESolidElement : public FEElement
{
public:
	//! default constructor
	FESolidElement(){}

	//! copy constructor
	FESolidElement(const FESolidElement& el);

	//! assignment operator
	FESolidElement& operator = (const FESolidElement& el);

	//! set the element traits
	void SetTraits(FEElementTraits* pt) override;

	double gr(int n) const { return ((FESolidElementTraits*)(m_pT))->gr[n]; }	// integration point coordinate r
	double gs(int n) const { return ((FESolidElementTraits*)(m_pT))->gs[n]; }	// integration point coordinate s
	double gt(int n) const { return ((FESolidElementTraits*)(m_pT))->gt[n]; }	// integration point coordinate t

	double* GaussWeights() const { return &((FESolidElementTraits*)(m_pT))->gw[0]; }			// weights of integration points

	double* Gr(int n) const { return ((FESolidElementTraits*)(m_pT))->Gr[n]; }	// shape function derivative to r
	double* Gs(int n) const { return ((FESolidElementTraits*)(m_pT))->Gs[n]; }	// shape function derivative to s
	double* Gt(int n) const { return ((FESolidElementTraits*)(m_pT))->Gt[n]; }	// shape function derivative to t

	double* Grr(int n) const { return ((FESolidElementTraits*)(m_pT))->Grr[n]; }	// shape function 2nd derivative to rr
	double* Gsr(int n) const { return ((FESolidElementTraits*)(m_pT))->Gsr[n]; }	// shape function 2nd derivative to sr
	double* Gtr(int n) const { return ((FESolidElementTraits*)(m_pT))->Gtr[n]; }	// shape function 2nd derivative to tr
	
	double* Grs(int n) const { return ((FESolidElementTraits*)(m_pT))->Grs[n]; }	// shape function 2nd derivative to rs
	double* Gss(int n) const { return ((FESolidElementTraits*)(m_pT))->Gss[n]; }	// shape function 2nd derivative to ss
	double* Gts(int n) const { return ((FESolidElementTraits*)(m_pT))->Gts[n]; }	// shape function 2nd derivative to ts
	
	double* Grt(int n) const { return ((FESolidElementTraits*)(m_pT))->Grt[n]; }	// shape function 2nd derivative to rt
	double* Gst(int n) const { return ((FESolidElementTraits*)(m_pT))->Gst[n]; }	// shape function 2nd derivative to st
	double* Gtt(int n) const { return ((FESolidElementTraits*)(m_pT))->Gtt[n]; }	// shape function 2nd derivative to tt

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t) const { ((FESolidElementTraits*)(m_pT))->shape_fnc(H, r, s, t); }

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) const { ((FESolidElementTraits*)(m_pT))->shape_deriv(Hr, Hs, Ht, r, s, t); }

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t) const { ((FESolidElementTraits*)(m_pT))->shape_deriv2(Hrr, Hss, Htt, Hrs, Hst, Hrt, r, s, t); }

	vec3d evaluate(vec3d* v, double r, double s, double t) const;

public:
	vector<bool>    m_bitfc;    //!< flag for interface nodes
	vector<mat3d>	m_J0i;		//!< inverse of reference Jacobian
};

//-----------------------------------------------------------------------------
//!  This class defines a surface element

class FECORE_API FESurfaceElement : public FEElement
{
public:
	FESurfaceElement();

	FESurfaceElement(const FESurfaceElement& el);

	FESurfaceElement& operator = (const FESurfaceElement& el);

	virtual void SetTraits(FEElementTraits* pt) override;

	double* GaussWeights() { return &((FESurfaceElementTraits*)(m_pT))->gw[0]; }			// weights of integration points
	double gr(int n) const { return ((FESurfaceElementTraits*)(m_pT))->gr[n]; }	// integration point coordinate r
	double gs(int n) const { return ((FESurfaceElementTraits*)(m_pT))->gs[n]; }	// integration point coordinate  s

	double* Gr(int n) const { return ((FESurfaceElementTraits*)(m_pT))->Gr[n]; }	// shape function derivative to r
	double* Gs(int n) const { return ((FESurfaceElementTraits*)(m_pT))->Gs[n]; }	// shape function derivative to s

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

    //! return number of edges
    int facet_edges();
    
    //! return node list of edge
    void facet_edge(int j, int* en);

	//! serialize
	void Serialize(DumpStream& ar) override;
    
public:
    //! local ID of surface element
	int		m_lid;

	// indices of solid or shell element this surface is a face of
	// For solids, a surface element can be connected to two elements 
	// if the surface is an inside surface. For boundary surfaces
	// the second element index is -1. 
	FEElement*		m_elem[2];				
};

//-----------------------------------------------------------------------------
//!  This class defines the shell element. 

//! A shell element is similar to a surface
//! element except that it has a thickness. 

class FECORE_API FEShellElement : public FEElement
{
public:
	FEShellElement();

	//! copy constructor
	FEShellElement(const FEShellElement& el);

	//! assignment operator
	FEShellElement& operator = (const FEShellElement& el);

	virtual void SetTraits(FEElementTraits* ptraits) override;

    double gr(int n) { return ((FEShellElementTraits*)(m_pT))->gr[n]; }
    double gs(int n) { return ((FEShellElementTraits*)(m_pT))->gs[n]; }
    double gt(int n) { return ((FEShellElementTraits*)(m_pT))->gt[n]; }
    
	double* GaussWeights() { return &((FEShellElementTraits*)(m_pT))->gw[0]; }	// weights of integration points

	double* Hr(int n) { return ((FEShellElementTraits*)(m_pT))->Hr[n]; }	// shape function derivative to r
	double* Hs(int n) { return ((FEShellElementTraits*)(m_pT))->Hs[n]; }	// shape function derivative to s

    //! values of shape functions
    void shape_fnc(double* H, double r, double s) const { ((FEShellElementTraits*)(m_pT))->shape_fnc(H, r, s); }
    
    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double r, double s) const { ((FEShellElementTraits*)(m_pT))->shape_deriv(Hr, Hs, r, s); }
    
    //! serialize data associated with this element
    void Serialize(DumpStream &ar) override;

public:
	vector<double>	m_h0;	//!< initial shell thicknesses
	vector<double>	m_ht;	//!< current shell thickness

    // indices of solid elements this shell element is attached to.
    // the first element is attached to the back of the shell
    // and the second element is attached to the front.
    // the index is -1 if no solid is attached on that side.
    int        m_elem[2];
};

//-----------------------------------------------------------------------------
// Shell element used by old shell formulation
class FECORE_API FEShellElementOld : public FEShellElement
{
public:
	FEShellElementOld();

	//! copy constructor
	FEShellElementOld(const FEShellElementOld& el);

	//! assignment operator
	FEShellElementOld& operator = (const FEShellElementOld& el);

	// set the element traits class
	void SetTraits(FEElementTraits* ptraits) override;

	//! serialize data associated with this element
	void Serialize(DumpStream &ar) override;

public:
	vector<vec3d>	m_D0;	//!< initial shell directors
};

//-----------------------------------------------------------------------------
// Shell element used by new shell formulations
class FECORE_API FEShellElementNew : public FEShellElement
{
public:
	FEShellElementNew();

	//! copy constructor
	FEShellElementNew(const FEShellElementNew& el);

	//! assignment operator
	FEShellElementNew& operator = (const FEShellElementNew& el);

	// set the element traits class
	void SetTraits(FEElementTraits* ptraits) override;

	//! serialize data associated with this element
	void Serialize(DumpStream &ar) override;
	
public: // EAS parameters

	matrix          m_Kaai;
	matrix          m_fa;
	matrix          m_alpha;
	matrix          m_alphat;
	matrix          m_alphai;
	vector<matrix>  m_Kua;
	vector<matrix>  m_Kwa;
	vector<mat3ds>  m_E;
};

//-----------------------------------------------------------------------------

class FECORE_API FETrussElement : public FEElement
{
public:
	FETrussElement();

	FETrussElement(const FETrussElement& el);

	FETrussElement& operator = (const FETrussElement& el);

public:
	double	m_a0;	// cross-sectional area
};

//-----------------------------------------------------------------------------
//! Discrete element class

class FECORE_API FEDiscreteElement : public FEElement
{
public:
	FEDiscreteElement(){}
	FEDiscreteElement(const FEDiscreteElement& e);
	FEDiscreteElement& operator = (const FEDiscreteElement& e);
};

//-----------------------------------------------------------------------------
//!  This class defines a 2D element
class FECORE_API FEElement2D : public FEElement
{
public:
	//! default constructor
	FEElement2D(){}

	//! copy constructor
	FEElement2D(const FEElement2D& el);

	//! assignment operator
	FEElement2D& operator = (const FEElement2D& el);

	double* GaussWeights() { return &((FE2DElementTraits*)(m_pT))->gw[0]; }			// weights of integration points

	double* Hr(int n) { return ((FE2DElementTraits*)(m_pT))->Gr[n]; }	// shape function derivative to r
	double* Hs(int n) { return ((FE2DElementTraits*)(m_pT))->Gs[n]; }	// shape function derivative to s

    double* Hrr(int n) { return ((FE2DElementTraits*)(m_pT))->Grr[n]; }	// shape function 2nd derivative to rr
    double* Hsr(int n) { return ((FE2DElementTraits*)(m_pT))->Gsr[n]; }	// shape function 2nd derivative to sr
    
    double* Hrs(int n) { return ((FE2DElementTraits*)(m_pT))->Grs[n]; }	// shape function 2nd derivative to rs
    double* Hss(int n) { return ((FE2DElementTraits*)(m_pT))->Gss[n]; }	// shape function 2nd derivative to ss
    
	//! values of shape functions
	void shape_fnc(double* H, double r, double s) { ((FE2DElementTraits*)(m_pT))->shape(H, r, s); }

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) { ((FE2DElementTraits*)(m_pT))->shape_deriv(Hr, Hs, r, s); }
};

//-----------------------------------------------------------------------------
class FECORE_API FELineElement : public FEElement
{
public:
	FELineElement();

	FELineElement(const FELineElement& el);

	FELineElement& operator = (const FELineElement& el);

	void SetTraits(FEElementTraits* pt);

	void Serialize(DumpStream& ar);

public:
	int	m_lid;	//!< local ID
};

#endif // !defined(AFX_FEELEMENT_H__2EE38101_58E2_4FEB_B214_BB71B6FB15FB__INCLUDED_)
