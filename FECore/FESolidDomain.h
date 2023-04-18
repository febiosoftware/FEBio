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
#include "FEDomain.h"
#include "FEDofList.h"
#include "FELinearSystem.h"
#include "FESolidElement.h"

//-----------------------------------------------------------------------------
// This typedef defines a surface integrand. 
// It evaluates the function at surface material point mp, and returns the value
// it the val vector. The size of the vector is determined by the field variable
// that is being integrated and is already set when the integrand is called.
// This is used in the FESurface::LoadVector function.
typedef std::function<void(FEMaterialPoint& mp, int node_a, std::vector<double>& val)> FEVolumeVectorIntegrand;

typedef std::function<void(FEMaterialPoint& mp, int node_a, int node_b, matrix& val)> FEVolumeMatrixIntegrand;

//-----------------------------------------------------------------------------
//! abstract base class for 3D volumetric elements
class FECORE_API FESolidDomain : public FEDomain
{
    FECORE_SUPER_CLASS(FESOLIDDOMAIN_ID)
    FECORE_BASE_CLASS(FESolidDomain)

public:
    //! constructor
    FESolidDomain(FEModel* pfem);
    
    //! create storage for elements
	bool Create(int nsize, FE_Element_Spec espec) override;

    //! return nr of elements
	int Elements() const override;

	//! initialize element data
	bool Init() override;
    
    //! reset data (overridden from FEDomain)
    void Reset() override;
    
    //! copy data from another domain (overridden from FEDomain)
    void CopyFrom(FEMeshPartition* pd) override;

    //! element access
	FESolidElement& Element(int n);
    FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

    int GetElementType() const { return m_Elem[0].Type(); }
    
    int GetElementShape() const { return m_Elem[0].Shape(); }

	FE_Element_Spec GetElementSpec() const;
    
    //! find the element in which point y lies
    FESolidElement* FindElement(const vec3d& y, double r[3]);

	//! find the element in which point y lies (reference configuration)
	FESolidElement* FindReferenceElement(const vec3d& y, double r[3]);

	//! Project a point to an element and return natural coordinates
	void ProjectToElement(FESolidElement& el, const vec3d& p, double r[3]);

	//! Project a point to an element in the reference frame and return natural coordinates
	//! returns true if the point lies in the element
	bool ProjectToReferenceElement(FESolidElement& el, const vec3d& p, double r[3]);

    //! Calculate deformation gradient at integration point n
    double defgrad(FESolidElement& el, mat3d& F, int n);
	double defgrad(FESolidElement& el, mat3d& F, int n, vec3d* r);

    //! Calculate deformation gradient at integration point n
    double defgrad(FESolidElement& el, mat3d& F, double r, double s, double t);
    
    //! Calculate deformation gradient at integration point n at previous time
    double defgradp(FESolidElement& el, mat3d& F, int n);
    
    //! Calculate GradJ at integration point n at current time
    vec3d GradJ(FESolidElement& el, int n);
    
    //! Calculate GradJ at integration point n at prev time
    vec3d GradJp(FESolidElement& el, int n);
    
    //! Calculate GradJ at integration point n at current time
    tens3dls Gradb(FESolidElement& el, int n);
    
    //! Calculate GradJ at integration point n at prev time
    tens3dls Gradbp(FESolidElement& el, int n);
    
    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjac0(const FESolidElement& el, double J[3][3], int n);
    
    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjac0(const FESolidElement& el, double J[3][3], double r, double s, double t);
    
    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjac0(const FESolidElement& el, double r, double s, double t, mat3d& J);
    
    //! calculate inverse jacobian matrix w.r.t. current frame
    double invjact(FESolidElement& el, double J[3][3], int n);
	double invjact(FESolidElement& el, double J[3][3], int n, const vec3d* r);

    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjact(FESolidElement& el, double J[3][3], double r, double s, double t);
    
    //! calculate inverse jacobian matrix w.r.t. intermediate frame
    double invjact(FESolidElement& el, double J[3][3], int n, const double alpha);
    
    //! calculate inverse jacobian matrix w.r.t. previous time
    double invjactp(FESolidElement& el, double J[3][3], int n);
    
    //! calculate gradient of function at integration points
    vec3d gradient(FESolidElement& el, double* fn, int n);
	vec3d gradient(FESolidElement& el, int order, double* fn, int n);

    //! calculate gradient of function at integration points
    vec3d gradient(FESolidElement& el, vector<double>& fn, int n);
    
    //! calculate spatial gradient of vector function at integration points
    mat3d gradient(FESolidElement& el, vec3d* fn, int n);
    
    //! calculate gradient of function at integration points at intermediate time
    vec3d gradient(FESolidElement& el, vector<double>& fn, int n, const double alpha);
    
    //! calculate spatial gradient of vector function at integration points at intermediate time
    mat3d gradient(FESolidElement& el, vec3d* fn, int n, const double alpha);
    
    //! calculate spatial gradient of vector function at integration points
    //! at previous time
    mat3d gradientp(FESolidElement& el, vec3d* fn, int n);
    
    //! calculate spatial gradient of vector function at integration points
    tens3dls gradient(FESolidElement& el, mat3ds* fn, int n);
    
    //! calculate spatial gradient of vector function at integration points
    //! at previous time
    tens3dls gradientp(FESolidElement& el, mat3ds* fn, int n);
    
    //! calculate material gradient of scalar function at integration points
    vec3d Gradient(FESolidElement& el, double* fn, int n);
    
    //! calculate material gradient of scalar function at integration points
    vec3d Gradient(FESolidElement& el, vector<double>& fn, int n);
    
    //! calculate material gradient of vector function at integration points
    mat3d Gradient(FESolidElement& el, vec3d* fn, int n);
    
    //! calculate material gradient of vector function at integration points
    tens3dls Gradient(FESolidElement& el, mat3ds* fn, int n);
    
    //! calculate jacobian in reference frame
    double detJ0(FESolidElement& el, int n);
    
    //! calculate jacobian in current frame
    double detJt(FESolidElement& el, int n);
    
    //! calculate jacobian in current frame
    double detJt(FESolidElement& el, int n, const double alpha);
    
    //! calculates covariant basis vectors in reference configuration at an integration point
    void CoBaseVectors0(FESolidElement& el, int j, vec3d g[3]);
    
    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FESolidElement& el, int j, vec3d g[3]);
    
    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FESolidElement& el, int j, vec3d g[3], const double alpha);
    
    //! calculates contravariant basis vectors in reference configuration at an integration point
    void ContraBaseVectors0(FESolidElement& el, int j, vec3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FESolidElement& el, int j, vec3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FESolidElement& el, int j, vec3d g[3], const double alpha);
    
    //! calculates parametric derivatives of covariant basis vectors at an integration point
    void CoBaseVectorDerivatives(FESolidElement& el, int j, vec3d dg[3][3]);
    
    //! calculates parametric derivatives of covariant basis vectors at an integration point
    void CoBaseVectorDerivatives0(FESolidElement& el, int j, vec3d dg[3][3]);
    
    //! calculates parametric derivatives of covariant basis vectors at an integration point
    void CoBaseVectorDerivatives(FESolidElement& el, int j, vec3d dg[3][3], const double alpha);
    
    //! calculates parametric derivatives of contravariant basis vectors at an integration point
    void ContraBaseVectorDerivatives(FESolidElement& el, int j, vec3d dg[3][3]);
    
    //! calculates parametric derivatives of contravariant basis vectors at an integration point
    void ContraBaseVectorDerivatives0(FESolidElement& el, int j, vec3d dg[3][3]);
    
    //! calculates parametric derivatives of contravariant basis vectors at an integration point
    void ContraBaseVectorDerivatives(FESolidElement& el, int j, vec3d dg[3][3], const double alpha);
    
    //! calculate the laplacian of a vector function at an integration point
    vec3d lapvec(FESolidElement& el, vec3d* fn, int n);
    
    //! calculate the gradient of the divergence of a vector function at an integration point
    vec3d gradivec(FESolidElement& el, vec3d* fn, int n);
    
    //! calculate the transpose of the gradient of the shape function gradients at an integration point
    void gradTgradShape(FESolidElement& el, int j, vector<mat3d>& mn);
    
    //! calculate spatial gradient of shapefunctions at integration point (returns Jacobian determinant)
    double ShapeGradient(FESolidElement& el, int n, vec3d* GradH);
    
    //! calculate spatial gradient of shapefunctions at integration point (returns Jacobian determinant)
    double ShapeGradient(FESolidElement& el, int n, vec3d* GradH, const double alpha);
    
    //! calculate spatial gradient of shapefunctions at integration point in reference frame (returns Jacobian determinant)
    double ShapeGradient0(FESolidElement& el, int n, vec3d* GradH);
    
    //! calculate spatial gradient of shapefunctions at integration point (returns Jacobian determinant)
    double ShapeGradient(FESolidElement& el, double r, double s, double t, vec3d* GradH);
    
    //! calculate spatial gradient of shapefunctions at integration point in reference frame (returns Jacobian determinant)
    double ShapeGradient0(FESolidElement& el, double r, double s, double t, vec3d* GradH);

	//! calculate the volume of an element in reference frame
	double Volume(FESolidElement& el);

	//! calculate the volume of an element in current frame
	double CurrentVolume(FESolidElement& el);

public:
	//! get the current nodal coordinates
	void GetCurrentNodalCoordinates(const FESolidElement& el, vec3d* rt);
	void GetCurrentNodalCoordinates(const FESolidElement& el, vec3d* rt, double alpha);

	//! get the reference nodal coordinates
	void GetReferenceNodalCoordinates(const FESolidElement& el, vec3d* r0);

	//! get the nodal coordinates at previous state
	void GetPreviousNodalCoordinates(const FESolidElement& el, vec3d* rp);

public:
	//! loop over elements
	void ForEachSolidElement(std::function<void(FESolidElement& el)> f);

	//! return the degrees of freedom of an element for this domain
	virtual int GetElementDofs(FESolidElement& el);

public:
	// Evaluate an integral over the domain and assemble into global load vector
	virtual void LoadVector(
		FEGlobalVector& R,			// the global vector to assembe the load vector in
		const FEDofList& dofList,	// the degree of freedom list
		FEVolumeVectorIntegrand f	// the actual integrand function
	);

	//! Evaluate the stiffness matrix of a load
	virtual void LoadStiffness(
		FELinearSystem& LS,			// The solver does the assembling
		const FEDofList& dofList_a,	// The degree of freedom list of node a
		const FEDofList& dofList_b,	// The degree of freedom list of node b
		FEVolumeMatrixIntegrand f	// the matrix function to evaluate
	);

protected:
    vector<FESolidElement>	m_Elem;		//!< array of elements
	FE_Element_Spec			m_elemSpec;	//!< the element spec

	FEDofList	m_dofU;
	FEDofList	m_dofSU;
};
