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



#pragma once
#include "matrix.h"
#include "vec3d.h"
#include "mat3d.h"
#include "fecore_enum.h"
#include <vector>

//-----------------------------------------------------------------------------
// Forward declaration of the FEElement class
class FEElement;
class FESolidElementShape;
class FESurfaceElementShape;

//-----------------------------------------------------------------------------
//! This class is the base class for all element trait's classes

class FECORE_API FEElementTraits
{
public:
	//! constructor
	FEElementTraits(int ni, int ne, FE_Element_Class c, FE_Element_Shape s, FE_Element_Type t);

	//! destructor
	virtual ~FEElementTraits(){}

	//! return the element class
	FE_Element_Class Class() const { return m_spec.eclass; }

	//! return the element shape
	FE_Element_Shape Shape() const { return m_spec.eshape; }

	//! return the element type
	FE_Element_Type Type() const { return m_spec.etype; }

	// project integration point data to nodes
	virtual void project_to_nodes(double* ai, double* ao) const {}
	virtual void project_to_nodes(vec3d*  ai, vec3d*  ao) const;
	virtual void project_to_nodes(mat3ds* ai, mat3ds* ao) const;
	virtual void project_to_nodes(mat3d*  ai, mat3d*  ao) const;

	virtual int ShapeFunctions(int order) { return m_neln; }

	int Faces() const { return m_faces; }

public:
	int m_nint;	//!< number of integration points
	int	m_neln;	//!< number of element nodes

	matrix m_H;	//!< shape function values at gausspoints.
				//!< The first index refers to the gauss-point,
				//!< the second index to the shape function

	std::vector<matrix> m_Hp;	//!< shape function values at gausspoints.
							//!< The first index refers to the gauss-point,
							//!< the second index to the shape function

	FE_Element_Spec	m_spec;	//!< element specs

protected:
	// number of faces of element
	int	m_faces;

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
class FECORE_API FESolidElementTraits : public FEElementTraits
{
public:
	//! constructor
	FESolidElementTraits(int ni, int ne, FE_Element_Shape es, FE_Element_Type et);

	//! initialize element traits data
	void init() override;

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);

	int ShapeFunctions(int order) override;

public:
	// gauss-point coordinates and weights
	std::vector<double> gr;
	std::vector<double> gs;
	std::vector<double> gt;
	std::vector<double> gw;

	// element shape class
	FESolidElementShape*				m_shape;
	std::vector<FESolidElementShape*>	m_shapeP; // shape classes for different order (some orders can by null)

	// local derivatives of shape functions at gauss points
	matrix m_Gr, m_Gs, m_Gt;
	std::vector<matrix>	m_Gr_p;
	std::vector<matrix>	m_Gs_p;
	std::vector<matrix>	m_Gt_p;

	// local second derivatives of shape functions at gauss points
	matrix Grr, Gsr, Gtr, Grs, Gss, Gts, Grt, Gst, Gtt;
};

//=============================================================================
class FECORE_API FESRISolidElementTraits
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

	//! initialize element traits data
	void init();

public:
	FEHex8_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_HEX8, et) {}
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

	void project_to_nodes(double* ai, double* ao) const override;

protected:
	matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
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

	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
// 8-node hexahedral element with uniform deformation gradient

class FEHex8G1 : public FEHex8_
{
public:
	enum { NINT = 1 };

public:
	FEHex8G1();

	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//! Base class for 4-node linear tetrahedrons
class FETet4_ : public FESolidElementTraits
{
public:
	enum { NELN = 4 };

	//! initialize element traits data
	void init();

public:
	FETet4_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_TET4, et) {}
};

//=============================================================================
// single Gauss point integrated tet element
class FETet4G1 : public FETet4_
{
public:
	enum { NINT = 1};

public:
	FETet4G1();

	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
// 4-node tetrahedral element using a 4-node Gaussian integration rule
class FETet4G4 : public FETet4_
{
public:
	enum { NINT = 4 };

public:
	FETet4G4();

	void project_to_nodes(double* ai, double* ao) const override;

protected:
	matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//! Base class for 5-node linear tetrahedrons
class FETet5_ : public FESolidElementTraits
{
public:
	enum { NELN = 5 };

	//! initialize element traits data
	void init();

public:
	FETet5_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_TET5, et) {}
};

//=============================================================================
// 5-node tetrahedral element using a 4-node Gaussian integration rule
class FETet5G4 : public FETet5_
{
public:
	enum { NINT = 4 };

public:
	FETet5G4();

	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//
//   FEPenta6
//
//=============================================================================

//=============================================================================
//! Base class for 6-node pentahedral "wedge" elements
class FEPenta6_ : public FESolidElementTraits
{
public:
	enum { NELN = 6 };

	//! initialize element traits data
	void init();

public:
	FEPenta6_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_PENTA6, et){}
};

//=============================================================================
// 6-node pentahedral elements with 6-point gaussian quadrature 
class FEPenta6G6 : public FEPenta6_
{
public:
	enum { NINT = 6 };

public:
	FEPenta6G6();

	void project_to_nodes(double* ai, double* ao) const override;

protected:
	matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
};

//=============================================================================
//
//   FEPenta15
//
//=============================================================================

//=============================================================================
//! Base class for 15-node quadratic pentahedral "wedge" elements
class FEPenta15_ : public FESolidElementTraits
{
public:
    enum { NELN = 15 };
    
public:
    FEPenta15_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_PENTA15, et){}
};

//=============================================================================
// 15-node pentahedral elements with 8-point gaussian quadrature
class FEPenta15G8 : public FEPenta15_
{
public:
    enum { NINT = 8 };
    
public:
    FEPenta15G8();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];
    
    matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
    matrix m_MT;
};

//=============================================================================
// 15-node pentahedral elements with 21-point gaussian quadrature
class FEPenta15G21 : public FEPenta15_
{
public:
    enum { NINT = 21 };
    
public:
    FEPenta15G21();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
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

	//! initialize element traits data
	void init();

public:
	FETet10_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_TET10, et){}
};

//=============================================================================
// 10-node tetrahedral element using a 4-node Gaussian integration rule
class FETet10G1 : public FETet10_
{
public:
	enum { NINT = 1 };

public:
	FETet10G1();

	void project_to_nodes(double* ai, double* ao) const override;

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

	void project_to_nodes(double* ai, double* ao) const override;

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

	void project_to_nodes(double* ai, double* ao) const override;

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

	void project_to_nodes(double* ai, double* ao) const override;
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
};

//=============================================================================
// 15-node tetrahedral element using a 8-node Gaussian integration rule
class FETet15G4 : public FETet15_
{
public:
	enum { NINT = 4 };

public:
	FETet15G4();

	void project_to_nodes(double* ai, double* ao) const override;

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

	void project_to_nodes(double* ai, double* ao) const override;

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

	void project_to_nodes(double* ai, double* ao) const override;

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

	void project_to_nodes(double* ai, double* ao) const override;

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
//   FETet20
//   
//=============================================================================

//=============================================================================
//! Base class for 20-node cubic tetrahedral elements
class FETet20_ : public FESolidElementTraits
{
public:
	enum { NELN = 20 };

public:
	FETet20_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_TET20, et){}
};

//=============================================================================
// 20-node tetrahedral element using a 15-point Gaussian integration rule
class FETet20G15 : public FETet20_
{
public:
	enum { NINT = 15 };

public:
	FETet20G15();

	void project_to_nodes(double* ai, double* ao) const override;

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
	FEHex20_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_HEX20, et){}

	//! initialize element traits data
	void init();

	int Nodes(int order);
};

//=============================================================================
// 20-node hexahedral element using a 2x2x2 Gaussian integration rule
class FEHex20G8 : public FEHex20_
{
public:
    enum { NINT = 8 };
    
public:
    FEHex20G8();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];
    
    matrix Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
    matrix MT;
};

//=============================================================================
// 20-node hexahedral element using a 3x3x3 Gaussian integration rule
class FEHex20G27 : public FEHex20_
{
public:
	enum { NINT = 27 };

public:
	FEHex20G27();

	void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    matrix Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
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
    
	//! initialize element traits data
	void init();

protected:
    matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 27-node hexahedral element using a 3x3x3 Gaussian integration rule
class FEHex27G27 : public FEHex27_
{
public:
	enum { NINT = 27 };

public:
	FEHex27G27();

	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
// 
// FEPyra5
//
//=============================================================================

//=============================================================================
//! Base class for 5-node pyramid element
class FEPyra5_ : public FESolidElementTraits
{
public:
	enum { NELN = 5 };

public:
	FEPyra5_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_PYRA5, et){}
};

//=============================================================================
// 5-node pyramid element using a 2x2x2 Gaussian integration rule
class FEPyra5G8: public FEPyra5_
{
public:
	enum { NINT = 8 };

public:
	FEPyra5G8();

	void project_to_nodes(double* ai, double* ao) const override;

protected:
	matrix	m_Ai;
};

//=============================================================================
//
// FEPyra13
//
//=============================================================================

//=============================================================================
//! Base class for 13-node pyramid element
class FEPyra13_ : public FESolidElementTraits
{
public:
    enum { NELN = 13 };
    
public:
    FEPyra13_(int ni, FE_Element_Type et) : FESolidElementTraits(ni, NELN, ET_PYRA13, et){}
};

//=============================================================================
// 13-node pyramid element using a 2x2x2 Gaussian integration rule
class FEPyra13G8: public FEPyra13_
{
public:
    enum { NINT = 8 };
    
public:
    FEPyra13G8();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    matrix    m_Ai;
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
class FECORE_API FESurfaceElementTraits : public FEElementTraits
{
public:
	FESurfaceElementTraits(int ni, int ne, FE_Element_Shape es, FE_Element_Type et);

	// initialization
	void init();

	// shape functions at (r,s)
	void shape_fnc(double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(double* Gr, double* Gs, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s);

	// shape functions at (r,s)
	void shape_fnc(int order, double* H, double r, double s);

	// shape function derivatives at (r,s)
	void shape_deriv(int order, double* Gr, double* Gs, double r, double s);

public:
	// gauss-point coordinates and weights
	std::vector<double> gr;
	std::vector<double> gs;
	std::vector<double> gw;

	// element shape class
	FESurfaceElementShape*				m_shape;
	std::vector<FESurfaceElementShape*>	m_shapeP; // shape classes for different order (some orders can by null)

	// local derivatives of shape functions at gauss points
	matrix Gr, Gs;

	// local derivatives of shape functions at gauss points, for different interpolation order
	std::vector<matrix>	Gr_p;
	std::vector<matrix>	Gs_p;
    
    // parametric coordinates of element center
    double cr;
    double cs;
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

	void init() override;

public:
	//! constructor
	FEQuad4_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, ET_QUAD4, et){}
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
	void project_to_nodes(double* ai, double* ao) const override;

protected:
	matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
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
	void project_to_nodes(double* ai, double* ao) const override;
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

	// initialization 
	void init() override;
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
	void project_to_nodes(double* ai, double* ao) const override;
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
	void project_to_nodes(double* ai, double* ao) const override;

protected:
	matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
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
	void project_to_nodes(double* ai, double* ao) const override;

protected:
	matrix	m_Ai;
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
	void project_to_nodes(double* ai, double* ao) const override;
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

	// initialization 
	void init() override;
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
	void project_to_nodes(double* ai, double* ao) const override;
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
	void project_to_nodes(double* ai, double* ao) const override;
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
	void project_to_nodes(double* ai, double* ao) const override;

private:
	matrix	m_Ai;
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
	void project_to_nodes(double* ai, double* ao) const override;
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
	void project_to_nodes(double* ai, double* ao) const override;
};

//=============================================================================
//
//   FETri6m
//   
//=============================================================================

//=============================================================================
// Base class for 6-noded quadratic triangles with modified shape functions
/*
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
	void project_to_nodes(double* ai, double* ao) const override;

private:
	matrix	m_Ai;
};
*/

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

	// initialization
	void init() override;
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
	void project_to_nodes(double* ai, double* ao) const override;

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
	void project_to_nodes(double* ai, double* ao) const override;

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
	void project_to_nodes(double* ai, double* ao) const override;

private:
	matrix	m_Ai;
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
	void project_to_nodes(double* ai, double* ao) const override;

private:
	matrix	Ai;
};

//=============================================================================
//
//   FETri10
//   
//=============================================================================

//=============================================================================
// Base class for 10-noded cubic triangles
class FETri10_ : public FESurfaceElementTraits
{
public:
	enum { NELN = 10 };

public:
	FETri10_(int ni, FE_Element_Type et) : FESurfaceElementTraits(ni, NELN, ET_TRI10, et){}

	// initialization
	void init() override;
};

//=============================================================================
//  10-node triangular element with 7-point gaussian quadrature
//
class FETri10G7 : public FETri10_
{
public:
	enum { NINT = 7 };

public:
	// constructor
	FETri10G7();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	matrix	m_Ai;
};


//=============================================================================
//  10-node triangular element with 12-point gaussian quadrature
//
class FETri10G12 : public FETri10_
{
public:
	enum { NINT = 12 };

public:
	// constructor
	FETri10G12();

	// project integration point data to nodes
	void project_to_nodes(double* ai, double* ao) const override;

private:
	matrix	m_Ai;
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

	// initialization
	void init() override;
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
	void project_to_nodes(double* ai, double* ao) const override;

private:
	matrix	m_Ai;
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
    void project_to_nodes(double* ai, double* ao) const override;
    
private:
    matrix	m_Ai;
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

	// initialization
	void init() override;
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
	void project_to_nodes(double* ai, double* ao) const override;

private:
	matrix	m_Ai;
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
    void project_to_nodes(double* ai, double* ao) const override;
    
private:
    matrix	m_Ai;
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

    void init();
    
    //! values of shape functions
    virtual void shape_fnc(double* H, double r, double s) = 0;
    
    //! values of shape function derivatives
    virtual void shape_deriv(double* Hr, double* Hs, double r, double s) = 0;
    
public:
	// gauss-point coordinates and weights
	std::vector<double> gr;
	std::vector<double> gs;
	std::vector<double> gt;
	std::vector<double> gw;

	// local derivatives of shape functions at gauss points
	matrix Hr, Hs;
    
};

//=============================================================================
// 4-node quadrilateral elements
//
class FEShellQuad4_ : public FEShellElementTraits
{
public:
    enum { NELN = 4 };
    
public:
    FEShellQuad4_(int ni, FE_Element_Type et) : FEShellElementTraits(ni, NELN, ET_QUAD4, et) {}
    
public:
    //! values of shape functions
    void shape_fnc(double* H, double r, double s);
    
    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double r, double s);
    
};

//=============================================================================
// 4-node quadrilateral elements with 4*2-point gaussian quadrature
//
class FEShellQuad4G8 : public FEShellQuad4_
{
public:
    enum { NINT = 8 };
    
public:
    FEShellQuad4G8();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 4-node quadrilateral elements with 4*3-point gaussian quadrature
//
class FEShellQuad4G12 : public FEShellQuad4_
{
public:
    enum { NINT = 12 };
    
public:
    FEShellQuad4G12();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 3-node triangular elements
//
class FEShellTri3_ : public FEShellElementTraits
{
public:
    enum { NELN = 3 };
    
public:
    FEShellTri3_(int ni, FE_Element_Type et) : FEShellElementTraits(ni, NELN, ET_TRI3, et) {}
    
public:
    //! values of shape functions
    void shape_fnc(double* H, double r, double s);
    
    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double r, double s);
    
};

//=============================================================================
// 3-node triangular elements with 3*2-point gaussian quadrature
//
class FEShellTri3G6 : public FEShellTri3_
{
public:
    enum { NINT = 6 };
    
public:
    FEShellTri3G6();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 3-node triangular elements with 3*3-point gaussian quadrature
//
class FEShellTri3G9 : public FEShellTri3_
{
public:
    enum { NINT = 9 };
    
public:
    FEShellTri3G9();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 8-node quadrilateral elements
//
class FEShellQuad8_ : public FEShellElementTraits
{
public:
    enum { NELN = 8 };
    
public:
    FEShellQuad8_(int ni, FE_Element_Type et) : FEShellElementTraits(ni, NELN, ET_QUAD8, et) {}
    
public:
    //! values of shape functions
    void shape_fnc(double* H, double r, double s);
    
    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double r, double s);
    
};

//=============================================================================
// 8-node quadrilateral elements with 9*2-point gaussian quadrature
//
class FEShellQuad8G18 : public FEShellQuad8_
{
public:
    enum { NINT = 18 };
    
public:
    FEShellQuad8G18();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 8-node quadrilateral elements with 9*3-point gaussian quadrature
//
class FEShellQuad8G27 : public FEShellQuad8_
{
public:
    enum { NINT = 27 };
    
public:
    FEShellQuad8G27();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 6-node triangular elements
//
class FEShellTri6_ : public FEShellElementTraits
{
public:
    enum { NELN = 6 };
    
public:
    FEShellTri6_(int ni, FE_Element_Type et) : FEShellElementTraits(ni, NELN, ET_TRI6, et) {}
    
public:
    //! values of shape functions
    void shape_fnc(double* H, double r, double s);
    
    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double r, double s);
    
};

//=============================================================================
// 6-node triangular elements with 7*2-point gaussian quadrature
//
class FEShellTri6G14 : public FEShellTri6_
{
public:
    enum { NINT = 14 };
    
public:
    FEShellTri6G14();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
};

//=============================================================================
// 6-node triangular elements with 7*3-point gaussian quadrature
//
class FEShellTri6G21 : public FEShellTri6_
{
public:
    enum { NINT = 21 };
    
public:
    FEShellTri6G21();
    
    void project_to_nodes(double* ai, double* ao) const override;
    
protected:
    // use these integration points to project to nodes
    static int ni[NELN];

    matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data
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

public:
	// gauss-point coordinates and weights
	std::vector<double> gr;
	std::vector<double> gs;
	std::vector<double> gw;

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
	void project_to_nodes(double* ai, double* ao) const override;
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
	void project_to_nodes(double* ai, double* ao) const override;
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
	void project_to_nodes(double* ai, double* ao) const override;

protected:
	matrix m_Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 
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
	void project_to_nodes(double* ai, double* ao) const override;

private:
	matrix	m_Ai;
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
	void project_to_nodes(double* ai, double* ao) const override;

private:
	matrix	m_Ai;
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

public:
	std::vector<double> gr;	//!< integration point coordinates
	std::vector<double> gw;	//!< integration point weights

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
	void project_to_nodes(double* ai, double* ao) const override;
};
