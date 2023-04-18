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
#include "FEElementLibrary.h"
#include "FEElementTraits.h"
#include "FEMaterialPoint.h"
#include "fecore_enum.h"
#include "FEException.h"

class FEMesh;
//-----------------------------------------------------------------------------
class FEElementTraits;
class FEMeshPartition;

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
	std::vector<FEMaterialPoint*>	m_data;
};

//-----------------------------------------------------------------------------
//! Base class for all element classes

//! From this class the different element classes are derived.

class FECORE_API FEElement
{
public:
	enum {MAX_NODES     = 27};	// max nr of nodes
	enum {MAX_INTPOINTS = 27};	// max nr of integration points

	// Status flags. 
	enum Status {
		ACTIVE = 0x01
	};

public:
	//! default constructor
	FEElement();

	//! destructor
	virtual ~FEElement() {}

	//! get the element ID
	int GetID() const;

	//! set the element ID
	void SetID(int n);

	//! Get the element's material ID
	int GetMatID() const;

	//! Set the element's material ID
	void SetMatID(int id);

	//Get the mesh partition that contains this element
	FEMeshPartition * GetMeshPartition() const { return m_part; }

	//Set the mesh partition that contains this element
	void SetMeshPartition(FEMeshPartition* part){ m_part = part; }

	//! Set the Local ID
	void SetLocalID(int lid) { m_lid = lid; }

	//! Get the local ID
	int GetLocalID() const { return m_lid; }

	//! clear material point data
	void ClearData();

public:
	//! Set the type of the element
	void SetType(int ntype) { FEElementLibrary::SetElementTraits(*this, ntype); }

	//! Set the traits of an element
	virtual void SetTraits(FEElementTraits* ptraits);

	//! Get the element traits
	FEElementTraits* GetTraits() { return m_pT; }

	//! return number of nodes
	int Nodes() const { return m_pT->m_neln; }

	//! return the element class
	int Class() const { return m_pT->Class(); }

	//! return the element shape
	int Shape() const { return m_pT->Shape(); }

	//! return the type of element
	int Type() const { return m_pT->Type(); }

	//! return number of integration points
	int GaussPoints() const { return m_pT->m_nint; }

	//! shape function values
	double* H(int n) { return m_pT->m_H[n]; }
	const double* H(int n) const { return m_pT->m_H[n]; }

	//! return number of faces
	int Faces() const { return m_pT->Faces(); }

	//! return the nodes of the face
	int GetFace(int nface, int* nodeList) const;

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
	double Evaluate(int order, double* fn, int n);

	//! evaluate scale field at integration point
	double Evaluate(std::vector<double>& fn, int n);

	//! evaluate vector field at integration point
	vec2d Evaluate(vec2d* vn, int n);

	//! evaluate vector field at integration point
	vec3d Evaluate(vec3d* vn, int n);

	// see if this element has the node n
    bool HasNode(int n) const;

    // see if this element has the list of nodes n
    int HasNodes(int* n, const int ns) const;
    
	// find local element index of node n
    int FindNode(int n) const;

	// project data to nodes
	void project_to_nodes(double* ai, double* ao) const { m_pT->project_to_nodes(ai, ao); }
	void project_to_nodes(vec3d*  ai, vec3d*  ao) const { m_pT->project_to_nodes(ai, ao); }
	void project_to_nodes(mat3ds* ai, mat3ds* ao) const { m_pT->project_to_nodes(ai, ao); }
	void project_to_nodes(mat3d*  ai, mat3d*  ao) const { m_pT->project_to_nodes(ai, ao); }

	// evaluate scalar field at integration point using specific interpolation order
	double Evaluate(double* fn, int order, int n);

	int ShapeFunctions(int order);
	double* H(int order, int n);

public:
	void setStatus(unsigned int n) { m_status = n; }
	unsigned int status() const { return m_status; }
	bool isActive() const { return (m_status & ACTIVE); }
	void setActive() { m_status |= ACTIVE; }
	void setInactive() { m_status &= ~ACTIVE; }

protected:
	int		m_nID;		//!< element ID
	int		m_lid;		//!< local ID
	int		m_mat;		//!< material index
	unsigned int	m_status;	//!< element status
	FEMeshPartition * m_part;	//!< parent mesh partition

public:
	std::vector<int>		m_node;		//!< connectivity

	// This array stores the local node numbers, that is the node numbers
	// into the node list of a domain.
	std::vector<int>		m_lnode;	//!< local connectivity

public: 
	// NOTE: Work in progress
	// Elements can now also have degrees of freedom, only currently just one.
	// Like with nodes, a degree of freedom needs an equation number and a value
	// The equation number is in m_lm and the value is in m_val
	int		m_lm;	//!< equation number of element degree of freedom
	double	m_val;	//!< solution value of element degree of freedom

protected:
	FEElementState		m_State;	//!< element state data
	FEElementTraits*	m_pT;		//!< pointer to element traits
};

//-----------------------------------------------------------------------------

class FECORE_API FETrussElement : public FEElement
{
public:
	FETrussElement();

	FETrussElement(const FETrussElement& el);

	FETrussElement& operator = (const FETrussElement& el);

	void Serialize(DumpStream& ar) override;

public:
	double	m_a0;	// cross-sectional area
	double	m_lam;	// current stretch ratio
	double	m_tau;	// Kirchoff stress
	double	m_L0;	// initial length
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
};
