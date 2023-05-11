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
#include "mat2d.h"
#include "vec2d.h"
#include "matrix.h"
#include "FEMeshPartition.h"
#include "FENodeSet.h"
#include "FEDofList.h"
#include "FESurfaceElement.h"
#include "FENode.h"

//-----------------------------------------------------------------------------
class FEMesh;
class FENodeSet;
class FEFacetSet;
class FELinearSystem;

//-----------------------------------------------------------------------------
class FECORE_API FESurfaceMaterialPoint : public FEMaterialPoint
{
public:
	vec3d	dxr, dxs;		// tangent vectors at material point

	// return the surface element
	FESurfaceElement* SurfaceElement() { return (FESurfaceElement*)m_elem; }

	void Serialize(DumpStream& ar) override
	{
		FEMaterialPoint::Serialize(ar);
		ar & dxr & dxs;
	}
};

// helper class for describing shape functions at dofs in integration routines
struct FECORE_API FESurfaceDofShape
{
	int			index;			// local dof index
	double		shape;			// shape functions
	double		shape_deriv_r;	// shape function r-derivative
	double		shape_deriv_s;	// shape function s-derivative
};

//-----------------------------------------------------------------------------
// This typedef defines a surface integrand. 
// It evaluates the function at surface material point mp, and returns the value
// it the val vector. The size of the vector is determined by the field variable
// that is being integrated and is already set when the integrand is called.
// This is used in the FESurface::LoadVector function.
typedef std::function<void(FESurfaceMaterialPoint& mp, const FESurfaceDofShape& node_a, std::vector<double>& val)> FESurfaceVectorIntegrand;

typedef std::function<void(FESurfaceMaterialPoint& mp, const FESurfaceDofShape& node_a, const FESurfaceDofShape& node_b, matrix& val)> FESurfaceMatrixIntegrand;

//-----------------------------------------------------------------------------
//! Surface mesh

//! This class implements the basic functionality for an FE surface.
//! More specialized surfaces are derived from this class

class FECORE_API FESurface : public FEMeshPartition
{
	FECORE_SUPER_CLASS(FESURFACE_ID)
	FECORE_BASE_CLASS(FESurface)

public:
	//! default constructor
	FESurface(FEModel* fem);

	//! destructor
	virtual ~FESurface();

	//! initialize surface data structure
	bool Init() override;
	void InitSurface();
    
	//! creates surface
	void Create(int nsize, int elemType = -1);

	//! Build a surface from a facet set
	void Create(const FEFacetSet& set);

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! unpack an LM vector from a dof list
	void UnpackLM(const FESurfaceElement& el, const FEDofList& dofList, vector<int>& lm);
	
	//! Extract a node set from this surface
	FENodeList GetNodeList();

	//! Get a list of bools that indicate whether the corresponding node is on the boundary
	// TODO: Move to MeshPartition
	void GetBoundaryFlags(std::vector<bool>& boundary) const;
    
    //! Set alpha parameter for intermediate time
    void SetAlpha(const double alpha) { m_alpha = alpha; }

public:

	//! return number of surface elements
	int Elements() const override { return (int)m_el.size(); }

	//! return an element of the surface
	FESurfaceElement& Element(int i) { return m_el[i]; }

	//! return an element of the surface
	const FESurfaceElement& Element(int i) const { return m_el[i]; }

	//! returns reference to element
	FEElement& ElementRef(int n) override { return m_el[n]; }
	const FEElement& ElementRef(int n) const override { return m_el[n]; }

	//! find the solid or shell element of a surface element
	FEElement* FindElement(FESurfaceElement& el);

    //! for interface surfaces, find the index of both solid elements
    //! on either side of the interface
    void FindElements(FESurfaceElement& el);

	//! loop over all elements
	void ForEachSurfaceElement(std::function<void(FESurfaceElement& el)> f);

public:
	// Create material point data for this surface
	virtual FEMaterialPoint* CreateMaterialPoint();

	// update surface data
	void Update(const FETimeInfo& tp) override;

public:

	//! Project a node onto a surface element
	vec3d ProjectToSurface(FESurfaceElement& el, vec3d x, double& r, double& s);

	//! check to see if a point is on element
	bool IsInsideElement(FESurfaceElement& el, double r, double s, double tol = 0);

	//! See if a ray intersects an element
	bool Intersect(FESurfaceElement& el, vec3d r, vec3d n, double rs[2], double& g, double eps);

	//! Invert the surface
	void Invert();

	//! Get the spatial position given natural coordinates
	vec3d Position(FESurfaceElement& el, double r, double s);

    //! Get the spatial position of an integration point
    vec3d Position(FESurfaceElement& el, int n);
    
	//! Get the nodal coordinates of an element
	void NodalCoordinates(FESurfaceElement& el, vec3d* re);
    
    //! Determine if a face on this surface is pointing away or into a specified element
    double FacePointing(FESurfaceElement& se, FEElement& el);
    

public:
	//! calculate the reference surface area of a surface element
	double FaceArea(FESurfaceElement& el);

	//! calculate the current surface area of a surface element
	double CurrentFaceArea(FESurfaceElement& el);

	//! return the max element size
	double MaxElementSize();

	//! calculate the metric tensor in the current configuration
	mat2d Metric(FESurfaceElement& el, double r, double s);

    //! calculate the metric tensor at an integration point
    mat2d Metric(const FESurfaceElement& el, int n) const;
    
    //! calculate the metric tensor at an integration point at previous time
    mat2d MetricP(FESurfaceElement& el, int n);
    
	//! calculate the metric tensor in the reference configuration
	mat2d Metric0(FESurfaceElement& el, double r, double s);

	//! calculate the surface normal
	vec3d SurfaceNormal(FESurfaceElement& el, double r, double s) const;

	//! calculate the surface normal at an integration point
	vec3d SurfaceNormal(const FESurfaceElement& el, int n) const;

    //! calculate the nodal normals
    void UpdateNodeNormals();
    
    //! return the nodal normals
    vec3d NodeNormal(const int inode) { return m_nn[inode]; }
    
	//! calculate the global position of a point on the surface
	vec3d Local2Global(FESurfaceElement& el, double r, double s);

	//! calculate the global position of an integration point
	vec3d Local2Global(FESurfaceElement& el, int n);

    //! calculate the global position of a point on the surface at previous time
    vec3d Local2GlobalP(FESurfaceElement& el, double r, double s);
    
    //! calculate the global position of an integration point at previous time
    vec3d Local2GlobalP(FESurfaceElement& el, int n);
    
	//! calculates the covariant base vectors of a surface at an integration point
	void CoBaseVectors(const FESurfaceElement& el, int j, vec3d t[2]) const;

	//! calculates the covariant base vectors of a surface
	void CoBaseVectors(FESurfaceElement& el, double r, double s, vec3d t[2]);

	//! calculates covariant base vectors of a surface
	void CoBaseVectors0(FESurfaceElement& el, double r, double s, vec3d t[2]);

    //! calculates the covariant base vectors of a surface at an integration point at previoust time step
    void CoBaseVectorsP(FESurfaceElement& el, int j, vec3d t[2]);
    
    //! calculates contravariant base vectors of a surface  at an integration point
    void ContraBaseVectors(const FESurfaceElement& el, int j, vec3d t[2]) const;
    
    //! calculates the contravariant base vectors of a surface at an integration point at previoust time step
    void ContraBaseVectorsP(FESurfaceElement& el, int j, vec3d t[2]);
    
    //! calculates the parametric derivatives of covariant basis of a surface  at an integration point
    void CoBaseVectorDerivatives(const FESurfaceElement& el, int j, vec3d dg[2][2]) const;
    
    //! calculates the the parametric derivatives of covariant basis of a surface at an integration point at previoust time step
    void CoBaseVectorDerivativesP(FESurfaceElement& el, int j, vec3d dg[2][2]);
    
	//! calculates contravariant base vectors of a surface
	void ContraBaseVectors(FESurfaceElement& el, double r, double s, vec3d t[2]);

	//! calculates contravariant base vectors of a surface
	void ContraBaseVectors0(FESurfaceElement& el, double r, double s, vec3d t[2]);

	//! Jacobian in reference configuration for integration point n
	double jac0(FESurfaceElement& el, int n);

	//! Jacobian in reference configuration for integration point n (and returns normal)
	double jac0(const FESurfaceElement& el, int n, vec3d& nu);

    //! Interface status
    void SetInterfaceStatus(const bool bitfc) { m_bitfc = bitfc; }
    bool GetInterfaceStatus() { return m_bitfc; }

	//! Get the facet set that created this surface
	FEFacetSet* GetFacetSet() { return m_surf; }

public:
	// Get nodal reference coordinates 
	void GetReferenceNodalCoordinates(FESurfaceElement& el, vec3d* r0);

	// Get current coordinates
	void GetNodalCoordinates(FESurfaceElement& el, vec3d* rt);

	// Get current coordinates at intermediate configuration
	void GetNodalCoordinates(FESurfaceElement& el, double alpha, vec3d* rt);

	// Get the shell bottom flag
	bool IsShellBottom() const { return m_bshellb; }

	// Set the shell bottom flag
	void SetShellBottom(bool b) { m_bshellb = b; }

public:
	// Evaluate field variables
	double Evaluate(FESurfaceMaterialPoint& mp, int dof);
    double Evaluate(int nface, int dof);

public:
	//! Evaluate a load vector. 
	virtual void LoadVector(
		FEGlobalVector& R,				// The global vector into which the loads are assembled
		const FEDofList& dofList,		// The degree of freedom list
		bool breference,				// integrate over reference (true) or current (false) configuration
		FESurfaceVectorIntegrand f);	// the function that evaluates the integrand

	//! Evaluate the stiffness matrix of a load
	virtual void LoadStiffness(
		FELinearSystem& LS,			// The linear system does the assembling
		const FEDofList& dofList_a,	// The degree of freedom list of node a
		const FEDofList& dofList_b,	// The degree of freedom list of node b
		FESurfaceMatrixIntegrand f	// the matrix function to evaluate
	);

public:
	void CreateMaterialPointData();
    
protected:
	FEFacetSet*					m_surf;		//!< the facet set from which this surface is built
	vector<FESurfaceElement>	m_el;		//!< surface elements
    vector<vec3d>               m_nn;       //!< node normals
    bool                        m_bitfc;    //!< interface status
    double                      m_alpha;    //!< intermediate time fraction
	bool						m_bshellb;	//!< true if this surface is the bottom of a shell domain
};
