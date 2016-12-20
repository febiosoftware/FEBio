// FEElementTraits.h: interface for the FEElementTraits class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEELEMENTTRAITS_H__5AE1C578_7EC7_4C11_AC98_EBCCFD68B00C__INCLUDED_)
#define AFX_FEELEMENTTRAITS_H__5AE1C578_7EC7_4C11_AC98_EBCCFD68B00C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "matrix.h"
#include "vec3d.h"
#include "mat3d.h"
#include "FE_enum.h"
#include <vector>

//-----------------------------------------------------------------------------
// Forward declaration of the FEElement class
class FEElement;

//-----------------------------------------------------------------------------
//! This class is the base class for all element trait's classes

class FEElementTraits
{
public:
	//! constructor
	FEElementTraits(int ni, int ne, FE_Element_Class c, FE_Element_Shape s, FE_Element_Type t);

	//! destructor
	virtual ~FEElementTraits(){}

	//! return the element class
	FE_Element_Class Class() const { return spec.eclass; }

	//! return the element shape
	FE_Element_Shape Shape() const { return spec.eshape; }

	//! return the element type
	FE_Element_Type Type() const { return spec.etype; }

public:
	int nint;	//!< number of integration points
	int	neln;	//!< number of element nodes

	matrix H;	//!< shape function values at gausspoints.
				//!< The first index refers to the gauss-point,
				//!< the second index to the shape function

	FE_Element_Spec	spec;	//!< element specs

private:

	//! function to allocate storage for integration point data
	virtual void init() = 0;
};

//=============================================================================
//      S O L I D   E L E M E N T 
//
// This section defines a set of solid element formulation used in 3D finite
// element models.
//=============================================================================

//=============================================================================
//! This class defines the specific traits for solid elements and serves as
//! a base class for specific solid element formulations
//
class FESolidElementTraits : public FEElementTraits
{
public:
	//! constructor
	FESolidElementTraits(int ni, int ne, FE_Element_Shape es, FE_Element_Type et);

	//! initialize element traits data
	void init();

	//! values of shape functions
	virtual void shape_fnc(double* H, double r, double s, double t) = 0;

	//! values of shape function derivatives
	virtual void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) = 0;

	//! values of shape function second derivatives
	virtual void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t) = 0;

	//! project integration point data to nodes
	virtual void project_to_nodes(double* ai, double* ao) = 0;

public:
	// gauss-point coordinates and weights
	vector<double> gr;
	vector<double> gs;
	vector<double> gt;
	vector<double> gw;

	// local derivatives of shape functions at gauss points
	matrix Gr, Gs, Gt;

	// local second derivatives of shape functions at gauss points
	matrix Grr, Gsr, Gtr, Grs, Gss, Gts, Grt, Gst, Gtt;
};

//=============================================================================
class FESRISolidElementTraits
{
public:
	FESRISolidElementTraits() : m_pTRI(0) {}
	~FESRISolidElementTraits() { if (m_pTRI) delete m_pTRI; }
	FESolidElementTraits*	m_pTRI;
};

//=============================================================================
//! Base class for 8-node hexahedral elements
class FEHex8_ : public FESolidElementTraits
{
public:
	enum { NELN = 8 };

public:
	FEHex8_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_HEX8, et) {}

public:
	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
// 8-node hexahedral elements with 8-point gaussian quadrature
//
class FEHex8G8 : public FEHex8_
{
public:
	enum { NINT = 8 };

public:
	FEHex8G8();

	void project_to_nodes(double* ai, double* ao);

protected:
	matrix Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
// 8-node hexahedral elements with 6-point reduced integration rule
//
class FEHex8RI : public FEHex8_
{
public:
	enum { NINT = 6 };

public:
	FEHex8RI();

	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
// 8-node hexahedral element with uniform deformation gradient

class FEHex8G1 : public FEHex8_
{
public:
	enum { NINT = 1 };

public:
	FEHex8G1();

	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//! Base class for 4-node linear tetrahedrons
class FETet4_ : public FESolidElementTraits
{
public:
	enum { NELN = 4 };

public:
	FETet4_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_TET4, et) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
// single Gauss point integrated tet element
class FETet4G1 : public FETet4_
{
public:
	enum { NINT = 1};

public:
	FETet4G1();

	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
// 4-node tetrahedral element using a 4-node Gaussian integration rule
class FETet4G4 : public FETet4_
{
public:
	enum { NINT = 4 };

public:
	FETet4G4();

	void project_to_nodes(double* ai, double* ao);

protected:
	matrix Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//! Base class for 6-node pentahedral "wedge" elements
class FEPenta6_ : public FESolidElementTraits
{
public:
	enum { NELN = 6 };

public:
	FEPenta6_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_PENTA6, et){}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
// 6-node pentahedral elements with 6-point gaussian quadrature 
class FEPenta6G6 : public FEPenta6_
{
public:
	enum { NINT = 6 };

public:
	FEPenta6G6();

	void project_to_nodes(double* ai, double* ao);

protected:
	matrix Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//
//   FETet10
//   
//=============================================================================

//=============================================================================
//! Base class for 10-node quadratic tetrahedral elements
class FETet10_ : public FESolidElementTraits
{
public:
	enum { NELN = 10 };

public:
	FETet10_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_TET10, et){}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
// 10-node tetrahedral element using a 4-node Gaussian integration rule
class FETet10G1 : public FETet10_
{
public:
	enum { NINT = 1 };

public:
	FETet10G1();

	void project_to_nodes(double* ai, double* ao);

private:
	matrix Ai;
};

//=============================================================================
// 10-node tetrahedral element using a 4-node Gaussian integration rule
class FETet10G4 : public FETet10_
{
public:
	enum { NINT = 4 };

public:
	FETet10G4();

	void project_to_nodes(double* ai, double* ao);

private:
	matrix Ai;
};

//=============================================================================
// 10-node tetrahedral element using a 8-node Gaussian integration rule
class FETet10G8 : public FETet10_
{
public:
	enum { NINT = 8 };

public:
	FETet10G8();

	void project_to_nodes(double* ai, double* ao);


private:
	matrix N;
	matrix Ai;
};

//=============================================================================
class FETet10G4RI1 : public FETet10G4, public FESRISolidElementTraits
{
public:
	FETet10G4RI1();
};

//=============================================================================
class FETet10G8RI4 : public FETet10G8, public FESRISolidElementTraits
{
public:
	FETet10G8RI4();
};

//=============================================================================
// 10-node tetrahedral element using a 11-node Gauss-Lobatto integration rule
class FETet10GL11 : public FETet10_
{
public:
	enum { NINT = 11 };

public:
	FETet10GL11();

	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//
//   FETet15
//   
//=============================================================================

//=============================================================================
//! Base class for 15-node quadratic tetrahedral elements
class FETet15_ : public FESolidElementTraits
{
public:
	enum { NELN = 15 };

public:
	FETet15_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_TET15, et){}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
// 15-node tetrahedral element using a 8-node Gaussian integration rule
class FETet15G4 : public FETet15_
{
public:
	enum { NINT = 4 };

public:
	FETet15G4();

	void project_to_nodes(double* ai, double* ao);

private:
	matrix N;
	matrix Ai;
};

//=============================================================================
// 15-node tetrahedral element using a 8-node Gaussian integration rule
class FETet15G8 : public FETet15_
{
public:
	enum { NINT = 8 };

public:
	FETet15G8();

	void project_to_nodes(double* ai, double* ao);

private:
	matrix N;
	matrix Ai;
};

//=============================================================================
// 15-node tetrahedral element using a 11-point Gaussian integration rule
class FETet15G11 : public FETet15_
{
public:
	enum { NINT = 11 };

public:
	FETet15G11();

	void project_to_nodes(double* ai, double* ao);

private:
	matrix N;
	matrix Ai;
};

//=============================================================================
// 15-node tetrahedral element using a 15-point Gaussian integration rule
class FETet15G15 : public FETet15_
{
public:
	enum { NINT = 15 };

public:
	FETet15G15();

	void project_to_nodes(double* ai, double* ao);

private:
	matrix N;
	matrix Ai;
};

//=============================================================================
class FETet15G15RI4 : public FETet15G15, public FESRISolidElementTraits
{
public:
	FETet15G15RI4();
};


//=============================================================================
//
//   FEHex20
//   
//=============================================================================


//=============================================================================
//! Base class for 20-node quadratic hexahedral element
class FEHex20_ : public FESolidElementTraits
{
public:
	enum { NELN = 20 };

public:
	FEHex20_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_HEX20, et){}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
// 20-node hexahedral element using a 3x3x3 Gaussian integration rule
class FEHex20G27 : public FEHex20_
{
public:
	enum { NINT = 27 };

public:
	FEHex20G27();

	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//
//   FEHex27
//   
//=============================================================================


//=============================================================================
//! Base class for 27-node quadratic hexahedral element
class FEHex27_ : public FESolidElementTraits
{
public:
	enum { NELN = 27 };

public:
	FEHex27_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_HEX27, et){}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
// 27-node hexahedral element using a 3x3x3 Gaussian integration rule
class FEHex27G27 : public FEHex27_
{
public:
	enum { NINT = 27 };

public:
	FEHex27G27();

	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//    S U R F A C E   E L E M E N T S
//
// This section defines a set of surface element formulations for use in 3D
// finite element models. For specific, these elements are used to define
// the surface of 3D volume models.
//=============================================================================

//=============================================================================
// This class defines the traits for surface elements and serves as a
// base class for the specific surface element formulations.
class FESurfaceElementTraits : public FEElementTraits
{
public:
	FESurfaceElementTraits(int ni, int ne, FE_Element_Shape es, FE_Element_Type et);

	// initialization
	void init();

	// shape functions at (r,s)
	virtual void shape(double* H, double r, double s) = 0;

	// shape function derivatives at (r,s)
	virtual void shape_deriv(double* Gr, double* Gs, double r, double s) = 0;

	// shape function derivatives at (r,s)
	virtual void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s) = 0;

	// project integration point data to nodes
	virtual void project_to_nodes(double* ai, double* ao) = 0;

public:
	// gauss-point coordinates and weights
	vector<double> gr;
	vector<double> gs;
	vector<double> gw;

	// local derivatives of shape functions at gauss points
	matrix Gr, Gs;
};

//=============================================================================
//
//   FEQuad4
//   
//=============================================================================

//=============================================================================
// Base class for 4-node bilinear quadrilaterals
//
class FEQuad4_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 4 };

public:
	//! constructor
	FEQuad4_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, ET_QUAD4, et){}

	//! shape functions at (r,s)
	void shape(double* H, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
// 4-node quadrilateral elements with 4-point gaussian quadrature 
class FEQuad4G4 : public FEQuad4_
{
public:
	enum { NINT = 4 };

public:
	FEQuad4G4();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

protected:
	matrix Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
// 4-node quadrilateral elements with nodal quadrature 
class FEQuad4NI : public FEQuad4_
{
public:
	enum { NINT = 4 };

public:
	//! constructor
	FEQuad4NI();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//
//   FETri3
//   
//=============================================================================

//=============================================================================
//! Base class for linear triangles
class FETri3_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 3 };

public:
	//! constructor
	FETri3_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, ET_TRI3, et){}

	//! shape function at (r,s)
	void shape(double* H, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//!  3-node triangular element with 1-point gaussian quadrature
class FETri3G1 : public FETri3_
{
public:
	enum { NINT = 1 };

public:
	//! constructor
	FETri3G1();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//!  3-node triangular element with 3-point gaussian quadrature
class FETri3G3 : public FETri3_
{
public:
	enum { NINT = 3 };

public:
	//! constructor
	FETri3G3();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

protected:
	matrix Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//!  3-node triangular element with 7-point gaussian quadrature
class FETri3G7 : public FETri3_
{
public:
	enum { NINT = 7 };

public:
	//! constructor
	FETri3G7();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

protected:
	matrix	Ai;
};

//=============================================================================
//  3-node triangular element with nodal quadrature
class FETri3NI : public FETri3_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	FETri3NI();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//
//   FETri6
//   
//=============================================================================

//=============================================================================
// Base class for 6-noded quadratic triangles
class FETri6_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 6 };

public:
	FETri6_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, ET_TRI6, et){}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//  6-node triangular element with 3-point gaussian quadrature
//
class FETri6G3 : public FETri6_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	FETri6G3();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//  6-node triangular element with 4-point gaussian quadrature
//
class FETri6G4 : public FETri6_
{
public:
	enum { NINT = 4 };

public:
	// constructor
	FETri6G4();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//  6-node triangular element with 7-point gaussian quadrature
//
class FETri6G7 : public FETri6_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri6G7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

private:
	matrix	Ai;
};

//=============================================================================
//  6-node triangular element with 7-point Gauss-Lobatto quadrature
//
class FETri6GL7 : public FETri6_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri6GL7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//  6-node triangular element with 6-point nodal quadrature
//
class FETri6NI : public FETri6_
{
public:
	enum { NINT = 6 };

public:
	// constructor
	FETri6NI();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//
//   FETri6m
//   
//=============================================================================

//=============================================================================
// Base class for 6-noded quadratic triangles with modified shape functions
class FETri6m_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 6 };

public:
	FETri6m_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, ET_TRI6, et){}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
// 6-node triangular element (with modified shape functions)
// with 7-point gaussian quadrature
//
class FETri6mG7 : public FETri6m_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri6mG7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

private:
	matrix	Ai;
};

//=============================================================================
//
//   FETri7
//   
//=============================================================================

//=============================================================================
// Base class for 7-noded quadratic triangles
class FETri7_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 7 };

public:
	FETri7_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, ET_TRI7, et){}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//  7-node triangular element with 3-point gaussian quadrature
//
class FETri7G3 : public FETri7_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	FETri7G3();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

private:
	matrix	Ai;
};

//=============================================================================
//  7-node triangular element with 4-point gaussian quadrature
//
class FETri7G4 : public FETri7_
{
public:
	enum { NINT = 4 };

public:
	// constructor
	FETri7G4();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

private:
	matrix	Ai;
};


//=============================================================================
//  7-node triangular element with 7-point gaussian quadrature
//
class FETri7G7 : public FETri7_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri7G7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

private:
	matrix	Ai;
};

//=============================================================================
//  7-node triangular element with 7-point Gauss-Lobatto quadrature
//
class FETri7GL7 : public FETri7_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri7GL7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

private:
	matrix	Ai;
};

//=============================================================================
//
//   FEQuad8
//   
//=============================================================================

//=============================================================================
//! Base class for 8-node quadratic quadrilaterals
//
class FEQuad8_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 8 };

public:
	FEQuad8_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, ET_QUAD8, et) {}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//! class implementing 8-node quad quadrilateral with 9 integration points
//
class FEQuad8G9 : public FEQuad8_
{
public:
	enum { NINT = 9 };

	// constructor
	FEQuad8G9();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

private:
	matrix	Ai;
};

//=============================================================================
//! class implementing 8-node quad quadrilateral with 8 nodal integration points
//
class FEQuad8NI : public FEQuad8_
{
public:
    enum { NINT = 8 };
    
    // constructor
    FEQuad8NI();
    
    // project integration point data to nodes
    void project_to_nodes(double* ai, double* ao);
    
private:
    matrix	Ai;
};

//=============================================================================
//
//   FEQuad9
//   
//=============================================================================

//=============================================================================
//! Base class for 9-node quadratic quadrilaterals
//
class FEQuad9_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 9 };

public:
	FEQuad9_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, ET_QUAD9, et) {}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//! class implementing 9-node quad quadrilateral with 9 integration points
//
class FEQuad9G9 : public FEQuad9_
{
public:
	enum { NINT = 9 };

	// constructor
	FEQuad9G9();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

private:
	matrix	Ai;
};

//=============================================================================
//! class implementing 9-node quad quadrilateral with 9 nodal integration points
//
class FEQuad9NI : public FEQuad9_
{
public:
    enum { NINT = 9 };
    
    // constructor
    FEQuad9NI();
    
    // project integration point data to nodes
    void project_to_nodes(double* ai, double* ao);
    
private:
    matrix	Ai;
};

//=============================================================================
//     S H E L L   E L E M E N T S
//
// This section defines several shell formulations for use in 3D finite element
// analysis. 
//=============================================================================

//=============================================================================
// This class defines the specific for shell elements and serves as a base class
// for specific shell formulations
//
class FEShellElementTraits : public FEElementTraits
{
public:
	FEShellElementTraits(int ni, int ne, FE_Element_Shape es, FE_Element_Type et);

public:
	// gauss-point coordinates and weights
	vector<double> gr;
	vector<double> gs;
	vector<double> gt;
	vector<double> gw;

	// local derivatives of shape functions at gauss points
	matrix Hr, Hs;
};

//=============================================================================
// 4-node quadrilateral elements with 4*3-point gaussian quadrature 
//
class FEShellQuadElementTraits : public FEShellElementTraits
{
public:
	enum { NINT = 12 };
	enum { NELN = 4 };

public:
	FEShellQuadElementTraits() : FEShellElementTraits(NINT, NELN, ET_QUAD4, FE_SHELL_QUAD) { init(); }

	void init();
};

//=============================================================================
// 3-node triangular elements with 3*3-point gaussian quadrature 
//
class FEShellTriElementTraits : public FEShellElementTraits
{
public:
	enum { NINT = 9 };
	enum { NELN = 3 };

public:
	FEShellTriElementTraits() : FEShellElementTraits(NINT, NELN, ET_TRI3, FE_SHELL_TRI) { init(); }

	void init();
};

//=============================================================================
// 8-node quadrilateral elements with 9*3-point gaussian quadrature
//
class FEShellQuad8ElementTraits : public FEShellElementTraits
{
public:
    enum { NINT = 27 };
    enum { NELN = 8 };
    
public:
    FEShellQuad8ElementTraits() : FEShellElementTraits(NINT, NELN, ET_QUAD8, FE_SHELL_QUAD8) { init(); }
    
    void init();
};

//=============================================================================
// 6-node triangular elements with 7*3-point gaussian quadrature
//
class FEShellTri6ElementTraits : public FEShellElementTraits
{
public:
    enum { NINT = 21 };
    enum { NELN = 6 };
    
public:
    FEShellTri6ElementTraits() : FEShellElementTraits(NINT, NELN, ET_TRI6, FE_SHELL_TRI6) { init(); }
    
    void init();
};

//=============================================================================
//          T R U S S    E L E M E N T S
//
// This section defines truss elements for 3D analysis
//=============================================================================

//=============================================================================
class FETrussElementTraits : public FEElementTraits
{
public:
	enum { NINT = 1 };
	enum { NELN = 2 };

public:
	FETrussElementTraits() : FEElementTraits(NINT, NELN, FE_ELEM_TRUSS, ET_TRUSS2, FE_TRUSS) { init(); }

	void init();
};

//=============================================================================
//          D I S C R E T E    E L E M E N T S
//
// This section defines discrete elements for 3D analysis
//=============================================================================

//=============================================================================
class FEDiscreteElementTraits : public FEElementTraits
{
public:
	enum { NINT = 1 };
	enum { NELN = 2 };

public:
	FEDiscreteElementTraits() : FEElementTraits(NINT, NELN, FE_ELEM_DISCRETE, ET_DISCRETE, FE_DISCRETE) { init(); }

	void init() {}
};

//=============================================================================
//                      2 D   E L E M E N T S
//
// This section defines a set of solid element formulation used in 3D finite
// element models.
//=============================================================================

//=============================================================================
// This class defines the traits for 2D elements and serves as a
// base class for the specific 2D element formulations.
class FE2DElementTraits : public FEElementTraits
{
public:
	FE2DElementTraits(int ni, int ne, FE_Element_Shape es, FE_Element_Type et);

	// initialization
	void init();

	// shape functions at (r,s)
	virtual void shape(double* H, double r, double s) = 0;

	// shape function derivatives at (r,s)
	virtual void shape_deriv(double* Gr, double* Gs, double r, double s) = 0;

	// shape function derivatives at (r,s)
	virtual void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s) = 0;

	// project integration point data to nodes
	virtual void project_to_nodes(double* ai, double* ao) = 0;

public:
	// gauss-point coordinates and weights
	vector<double> gr;
	vector<double> gs;
	vector<double> gw;

	// local derivatives of shape functions at gauss points
	matrix Gr, Gs;
    
    // local second derivatives of shape functions at gauss points
    matrix Grr, Gsr, Grs, Gss;
};

//=============================================================================
//
//   FE2DTri3
//   
//=============================================================================

//=============================================================================
//! Base class for linear triangles
class FE2DTri3_ : public FE2DElementTraits
{
public:
	enum { NELN = 3 };

public:
	//! constructor
	FE2DTri3_(int ni, FE_Element_Type et) : FE2DElementTraits(ni, NELN, ET_TRI3, et){}

	//! shape function at (r,s)
	void shape(double* H, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//!  3-node triangular element with 1-point gaussian quadrature
class FE2DTri3G1 : public FE2DTri3_
{
public:
	enum { NINT = 1 };

public:
	//! constructor
	FE2DTri3G1();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//
//   FE2DTri6
//   
//=============================================================================

//=============================================================================
// Base class for 6-noded quadratic triangles
class FE2DTri6_ : public FE2DElementTraits
{
public:
	enum { NELN = 6 };

public:
	FE2DTri6_(int ni, FE_Element_Type et) : FE2DElementTraits(ni, NELN, ET_TRI6, et){}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//  6-node triangular element with 3-point gaussian quadrature
//
class FE2DTri6G3 : public FE2DTri6_
{
public:
	enum { NINT = 3 };

public:
	// constructor
	FE2DTri6G3();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);
};

//=============================================================================
//
//   FE2DQuad4
//   
//=============================================================================

//=============================================================================
// Base class for 4-node bilinear quadrilaterals
//
class FE2DQuad4_ : public FE2DElementTraits
{
public:
	enum { NELN = 4 };

public:
	//! constructor
	FE2DQuad4_(int ni, FE_Element_Type et) : FE2DElementTraits(ni, NELN, ET_QUAD4, et){}

	//! shape functions at (r,s)
	void shape(double* H, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	//! shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
// 4-node quadrilateral elements with 4-point gaussian quadrature 
class FE2DQuad4G4 : public FE2DQuad4_
{
public:
	enum { NINT = 4 };

public:
	FE2DQuad4G4();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

protected:
	matrix Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//
//   FE2DQuad8
//   
//=============================================================================

//=============================================================================
//! Base class for 8-node quadratic quadrilaterals
//
class FE2DQuad8_ : public FE2DElementTraits
{
public:
	enum { NELN = 8 };

public:
	FE2DQuad8_(int ni, FE_Element_Type et) : FE2DElementTraits(ni, NELN, ET_QUAD8, et) {}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//! class implementing 8-node quad quadrilateral with 9 integration points
//
class FE2DQuad8G9 : public FE2DQuad8_
{
public:
	enum { NINT = 9 };

	// constructor
	FE2DQuad8G9();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

private:
	matrix	Ai;
};

//=============================================================================
//
//   FE2DQuad9
//   
//=============================================================================

//=============================================================================
//! Base class for 9-node quadratic quadrilaterals
//
class FE2DQuad9_ : public FE2DElementTraits
{
public:
	enum { NELN = 9 };

public:
	FE2DQuad9_(int ni, FE_Element_Type et) : FE2DElementTraits(ni, NELN, ET_QUAD9, et) {}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
};

//=============================================================================
//! class implementing 9-node quad quadrilateral with 9 integration points
//
class FE2DQuad9G9 : public FE2DQuad9_
{
public:
	enum { NINT = 9 };

	// constructor
	FE2DQuad9G9();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);

private:
	matrix	Ai;
};

//=============================================================================
//                      L I N E   E L E M E N T S
//
// This section defines a set of element formulations used to describe edges. 
//=============================================================================

//=============================================================================
class FELineElementTraits : public FEElementTraits
{
public:
	FELineElementTraits(int ni, int ne, FE_Element_Shape es, FE_Element_Type et);

	// initialization
	void init();

	// shape functions at r
	virtual void shape(double* H, double r) = 0;

	// shape function derivatives at (r)
	virtual void shape_deriv(double* Gr, double r) = 0;

	// shape function second derivatives at (r)
	virtual void shape_deriv2(double* Grr, double r) = 0;

	// project integration point data to nodes
	virtual void project_to_nodes(double* ai, double* ao) = 0;

public:
	vector<double> gr;	//!< integration point coordinates
	vector<double> gw;	//!< integration point weights

	// local derivatives of shape functions at gauss points
	matrix Gr;
    
    // local second derivatives of shape functions at gauss points
    matrix Grr;
};

//=============================================================================
//
//   FELine2_
//   
//=============================================================================

//=============================================================================
//! Base class for two-point lines
class FELine2_ : public FELineElementTraits
{
public:
	enum { NELN = 2 };

public:
	//! constructor
	FELine2_(int ni, FE_Element_Type et) : FELineElementTraits(ni, NELN, ET_LINE2, et){}

	//! shape function at (r)
	void shape(double* H, double r);

	//! shape function derivatives at (r)
	void shape_deriv(double* Gr, double r);

	//! shape function derivatives at (r)
	void shape_deriv2(double* Grr, double r);
};

//=============================================================================
class FELine2G1 : public FELine2_
{
public:
	enum { NINT = 1 };

public:
	//! constructor
	FELine2G1();

	//! project integration point data to nodes
	void project_to_nodes(double* ai, double* ao);
};

#endif // !defined(AFX_FEELEMENTTRAITS_H__5AE1C578_7EC7_4C11_AC98_EBCCFD68B00C__INCLUDED_)
