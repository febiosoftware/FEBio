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
// Max nr of nodal degrees of freedom

// At this point the nodal dofs are used as follows:
//
#define DOF_X			0		// x-displacement
#define DOF_Y			1		// y-displacement
#define DOF_Z			2		// z-displacement
#define DOF_U			3		// x-rotation
#define DOF_V			4		// y-rotation
#define DOF_W			5		// z-rotation
#define DOF_P			6		// fluid pressure
#define DOF_RU			7		// rigid x-rotation
#define DOF_RV			8		// rigid y-rotation
#define DOF_RW			9		// rigid z-rotation
#define DOF_T			10		// temperature
#define DOF_C			11		// solute concentration
//
// The rotational degrees of freedom are only used for rigid nodes and shells.
// The fluid pressure is only used for poroelastic problems.
// The rigid rotational degrees of freedom are only used for rigid nodes and only during the creation of the stiffenss matrix
// The temperature is only used during heat-conduction problems
// The solute concentration is only used in solute transport problems.

//-----------------------------------------------------------------------------
// Forward declaration of the FEElement class
class FEElement;

//-----------------------------------------------------------------------------
//! This class is the base class for all element trait's classes

class FEElementTraits
{
public:
	//! constructor
	FEElementTraits(int ni, int ne, FE_Element_Type et, FE_Element_Class ec) : m_ntype(et), m_nclass(ec)
	{
		neln = ne;
		nint = ni;
		
		H.resize(ni, ne);
	}

	//! destructor
	virtual ~FEElementTraits(){}

public:
	int nint;	//!< number of integration points
	int	neln;	//!< number of element nodes

	matrix H;	//!< shape function values at gausspoints.

				//!< The first index refers to the gauss-point,
				//!< the second index to the shape function

	int m_ntype;	//!< type of element
	int	m_nclass;	//!< element class

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
	FESolidElementTraits(int ni, int ne, FE_Element_Type et);

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
	
	// data used when unpacking
	// TODO: Are we still using this?
	vector<mat3d>	m_Jt;		// jacobian
	vector<mat3d>	m_Jti;		// inverse jacobian
	vector<double>	m_detJt;	// jacobian determinant

	vector<mat3d>	m_J0;		// jacobian
	vector<mat3d>	m_J0i;		// inverse jacobian
	vector<double>	m_detJ0;	// jacobian determinant
};

//=============================================================================
//! Base class for 8-node hexahedral elements
class FEHex8_ : public FESolidElementTraits
{
public:
	enum { NELN = 8 };

public:
	FEHex8_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, et) {}

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
	FETet4_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, et) {}

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
	FEPenta6_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, et){}

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
	FETet10_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, et){}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
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
	FETet15_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, et){}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
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
	FEHex20_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, et){}

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
	FEHex27_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, et){}

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
	FESurfaceElementTraits(int ni, int ne, FE_Element_Type et);

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
	FEQuad4_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, et){}

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
	FETri3_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, et){}

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
	FETri6_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, et){}

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
	FETri7_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, et){}

	// shape function at (r,s)
	void shape(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);
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
	FEQuad8_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, et) {}

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
	FEQuad9_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, et) {}

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
	FEShellElementTraits(int ni, int ne, FE_Element_Type et) : FEElementTraits(ni, ne, et, FE_ELEM_SHELL) 
	{
		gr.resize(ni);
		gs.resize(ni);
		gt.resize(ni);
		gw.resize(ni);

		Hr.resize(ni, ne);
		Hs.resize(ni, ne);

		D0.resize(ne);
		Dt.resize(ne);

		m_Jt.resize(ni);
		m_Jti.resize(ni);
		m_detJt.resize(ni);

		m_J0.resize(ni);
		m_J0i.resize(ni);
		m_detJ0.resize(ni);
	}

public:
	// gauss-point coordinates and weights
	vector<double> gr;
	vector<double> gs;
	vector<double> gt;
	vector<double> gw;

	// directors
	vector<vec3d>	D0;	//!< initial directors
	vector<vec3d>	Dt;	//!< current directors

	// local derivatives of shape functions at gauss points
	matrix Hr, Hs;

	// data used when unpacking
	vector<mat3d>	m_Jt;		// jacobian
	vector<mat3d>	m_Jti;		// inverse jacobian
	vector<double>	m_detJt;	// jacobian determinant

	vector<mat3d>	m_J0;		// jacobian
	vector<mat3d>	m_J0i;		// inverse jacobian
	vector<double>	m_detJ0;	// jacobian determinant
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
	FEShellQuadElementTraits() : FEShellElementTraits(NINT, NELN, FE_SHELL_QUAD) { init(); }

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
	FEShellTriElementTraits() : FEShellElementTraits(NINT, NELN, FE_SHELL_TRI) { init(); }

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
	FETrussElementTraits() : FEElementTraits(NINT, NELN, FE_TRUSS, FE_ELEM_TRUSS) { init(); }

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
	FEDiscreteElementTraits() : FEElementTraits(NINT, NELN, FE_DISCRETE, FE_ELEM_DISCRETE) { init(); }

	void init() {}
};

#endif // !defined(AFX_FEELEMENTTRAITS_H__5AE1C578_7EC7_4C11_AC98_EBCCFD68B00C__INCLUDED_)
