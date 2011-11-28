// FEElementTraits.h: interface for the FEElementTraits class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEELEMENTTRAITS_H__5AE1C578_7EC7_4C11_AC98_EBCCFD68B00C__INCLUDED_)
#define AFX_FEELEMENTTRAITS_H__5AE1C578_7EC7_4C11_AC98_EBCCFD68B00C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "vec3d.h"
#include "matrix.h"
#include "mat3d.h"
#include "FE_enum.h"
#include <vector>

//-----------------------------------------------------------------------------
// Max nr of nodal degrees of freedom

#define MAX_CDOFS	2
//#define MAX_NDOFS	11+MAX_CDOFS
#define MAX_NDOFS	13

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
// The rigid rotational degrees of freedom are only used for rigid nodes and only during the creation of the stiffness matrix
// The temperature is only used during heat-conduction problems
// The solute concentration is only used in solute transport problems.

//-----------------------------------------------------------------------------
// Forward declaration of the FEElement class
class FEElement;

//-----------------------------------------------------------------------------
//! This class is the base class for all element trait's classes
//
class FEElementTraits
{
public:
	//! constructor
	FEElementTraits(int ni, int ne)
	{
		neln = ne;
		nint = ni;
		
		H.Create(ni, ne);
	}

	//! destructor
	virtual ~FEElementTraits(){}

public:
	int nint;	//!< number of integration points
	int	neln;	//!< number of element nodes

	matrix H;	//!< shape function values at gausspoints.

				//!< The first index refers to the gauss-point,
				//!< the second index to the shape function

	matrix Hi;	//!< inverse of H; useful for projection integr. point data to nodal data 


	int m_ntype;	//!< type of element

private:

	//! function to set values of previous variables
	virtual void init() = 0;
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FESolidElementTraits
//  This class defines the specific traits for solid elements
//

class FESolidElementTraits : public FEElementTraits
{
public:
	FESolidElementTraits(int ni, int ne) : FEElementTraits(ni, ne) 
	{
		gr.resize(ni);
		gs.resize(ni);
		gt.resize(ni);
		gw.resize(ni);

		Gr.Create(ni, ne);
		Gs.Create(ni, ne);
		Gt.Create(ni, ne);

		Grr.Create(ni, ne);
		Gsr.Create(ni, ne);
		Gtr.Create(ni, ne);
		
		Grs.Create(ni, ne);
		Gss.Create(ni, ne);
		Gts.Create(ni, ne);
		
		Grt.Create(ni, ne);
		Gst.Create(ni, ne);
		Gtt.Create(ni, ne);
	}

public:
	// gauss-point coordinates and weights
	std::vector<double> gr;
	std::vector<double> gs;
	std::vector<double> gt;
	std::vector<double> gw;

	// local derivatives of shape functions at gauss points
	matrix Gr, Gs, Gt;

	// local second derivatives of shape functions at gauss points
	matrix Grr, Gsr, Gtr, Grs, Gss, Gts, Grt, Gst, Gtt;
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FEHexElementTraits
// 8-node hexahedral elements with 8-point gaussian quadrature
//

class FEHexElementTraits : public FESolidElementTraits
{
public:
	enum { NINT = 8 };
	enum { NELN = 8 };

public:
	FEHexElementTraits() : FESolidElementTraits(NINT, NELN) { m_ntype = FE_HEX; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FERIHexElementTraits
// 8-node hexahedral elements with 6-point reduced integration rule
//

class FERIHexElementTraits : public FESolidElementTraits
{
public:
	enum { NINT = 6 };
	enum { NELN = 8 };

public:
	FERIHexElementTraits() : FESolidElementTraits(NINT, NELN) { m_ntype = FE_RIHEX; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FEUDFHexElementTraits
// 8-node hexahedral element with uniform deformation gradient

class FEUDFHexElementTraits : public FESolidElementTraits
{
public:
	enum { NINT = 1 };
	enum { NELN = 8 };

public:
	FEUDFHexElementTraits() : FESolidElementTraits(NINT, NELN) { m_ntype = FE_UDGHEX; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FETetElementTraits
// 4-node tetrahedral elements with 4-point gaussian quadrature
//

class FETetElementTraits : public FESolidElementTraits
{
public:
	enum { NINT = 4 };
	enum { NELN = 4 };

public:
	FETetElementTraits() : FESolidElementTraits(NINT, NELN) { m_ntype = FE_TET; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
// single Gauss point integrated tet element

class FEG1TetElementTraits : public FESolidElementTraits
{
public:
	enum { NINT = 1};
	enum { NELN = 4};

public:
	FEG1TetElementTraits() : FESolidElementTraits(NINT, NELN) { m_ntype = FE_TETG1; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FEPentaElementTraits
// 6-node pentahedral elements with 6-point gaussian quadrature 
//

class FEPentaElementTraits : public FESolidElementTraits
{
public:
	enum { NINT = 6 };
	enum { NELN = 6 };

public:
	FEPentaElementTraits() : FESolidElementTraits(NINT, NELN) { m_ntype = FE_PENTA; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FESurfaceElementTraits
//  This class defines the specific traits for shell elements
//

class FESurfaceElementTraits : public FEElementTraits
{
public:
	FESurfaceElementTraits(int ni, int ne) : FEElementTraits(ni, ne) 
	{
		gr.resize(ni);
		gs.resize(ni);
		gw.resize(ni);

		Gr.Create(ni, ne);
		Gs.Create(ni, ne);
	}

public:
	// gauss-point coordinates and weights
	std::vector<double> gr;
	std::vector<double> gs;
	std::vector<double> gw;

	// local derivatives of shape functions at gauss points
	matrix Gr, Gs;
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FEQuadElementTraits
// 4-node quadrilateral elements with 4-point gaussian quadrature 
//

class FEQuadElementTraits : public FESurfaceElementTraits
{
public:
	enum { NINT = 4 };
	enum { NELN = 4 };

public:
	FEQuadElementTraits() : FESurfaceElementTraits(NINT, NELN) { m_ntype = FE_QUAD; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FENIQuadElementTraits
// 4-node quadrilateral elements with nodal quadrature 
//

class FENIQuadElementTraits : public FESurfaceElementTraits
{
public:
	enum { NINT = 4 };
	enum { NELN = 4 };

public:
	FENIQuadElementTraits() : FESurfaceElementTraits(NINT, NELN) { m_ntype = FE_NIQUAD; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FETriElementTraits
//  3-node triangular element with 3-point gaussian quadrature
//

class FETriElementTraits : public FESurfaceElementTraits
{
public:
	enum { NINT = 3 };
	enum { NELN = 3 };

public:
	FETriElementTraits() : FESurfaceElementTraits(NINT, NELN) { m_ntype = FE_TRI; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FENITriElementTraits
//  3-node triangular element with nodal quadrature
//

class FENITriElementTraits : public FESurfaceElementTraits
{
public:
	enum { NINT = 3 };
	enum { NELN = 3 };

public:
	FENITriElementTraits() : FESurfaceElementTraits(NINT, NELN) { m_ntype = FE_NITRI; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FEShellElementTraits
//  This class defines the specific traits for shell elements
//

class FEShellElementTraits : public FEElementTraits
{
public:
	FEShellElementTraits(int ni, int ne) : FEElementTraits(ni, ne) 
	{
		gr.resize(ni);
		gs.resize(ni);
		gt.resize(ni);
		gw.resize(ni);

		Hr.Create(ni, ne);
		Hs.Create(ni, ne);
	}

public:
	// gauss-point coordinates and weights
	std::vector<double> gr;
	std::vector<double> gs;
	std::vector<double> gt;
	std::vector<double> gw;

	// local derivatives of shape functions at gauss points
	matrix Hr, Hs;
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FEShellQuadElementTraits
// 4-node quadrilateral elements with 4*3-point gaussian quadrature 
//

class FEShellQuadElementTraits : public FEShellElementTraits
{
public:
	enum { NINT = 12 };
	enum { NELN = 4 };

public:
	FEShellQuadElementTraits() : FEShellElementTraits(NINT, NELN) { m_ntype = FE_SHELL_QUAD; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
// CLASS: FEShellTriElementTraits
// 3-node triangular elements with 3*3-point gaussian quadrature 
//

class FEShellTriElementTraits : public FEShellElementTraits
{
public:
	enum { NINT = 9 };
	enum { NELN = 3 };

public:
	FEShellTriElementTraits() : FEShellElementTraits(NINT, NELN) { m_ntype = FE_SHELL_TRI; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
class FETrussElementTraits : public FEElementTraits
{
public:
	enum { NINT = 1 };
	enum { NELN = 2 };

public:
	FETrussElementTraits() : FEElementTraits(NINT, NELN) { m_ntype = FE_TRUSS; init(); }

	void init();
};

///////////////////////////////////////////////////////////////////////////////
class FEDiscreteElementTraits : public FEElementTraits
{
public:
	enum { NINT = 1 };
	enum { NELN = 2 };

public:
	FEDiscreteElementTraits() : FEElementTraits(NINT, NELN) { m_ntype = FE_DISCRETE; init(); }

	void init() {}
};

#endif // !defined(AFX_FEELEMENTTRAITS_H__5AE1C578_7EC7_4C11_AC98_EBCCFD68B00C__INCLUDED_)
