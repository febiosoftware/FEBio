#pragma once
#include "FEDomain.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
//! abstract base class for 3D volumetric elements
class FESolidDomain : public FEDomain
{
public:
    //! constructor
    FESolidDomain(FEModel* pfem);
    
    //! create storage for elements
    void Create(int nsize, int elemType) override;
    
    //! return nr of elements
    int Elements() const override { return (int)m_Elem.size(); }

	//! initialize element data
	bool Init() override;
    
    //! reset data (overridden from FEDomain)
    void Reset() override;
    
    //! copy data from another domain (overridden from FEDomain)
    void CopyFrom(FEDomain* pd) override;
    
    //! element access
    FESolidElement& Element(int n) { return m_Elem[n]; }
    FEElement& ElementRef(int n) override { return m_Elem[n]; }
    
    int GetElementType() const { return m_Elem[0].Type(); }
    
    int GetElementShape() const { return m_Elem[0].Shape(); }
    
    //! find the element in which point y lies
    FESolidElement* FindElement(vec3d y, double r[3]);
    
    //! Calculate deformation gradient at integration point n
    double defgrad(FESolidElement& el, mat3d& F, int n);
    
    //! Calculate deformation gradient at integration point n
    double defgrad(FESolidElement& el, mat3d& F, double r, double s, double t);
    
    //! Calculate deformation gradient at integration point n at previous time
    double defgradp(FESolidElement& el, mat3d& F, int n);
    
    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjac0(const FESolidElement& el, double J[3][3], int n);
    
    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjac0(const FESolidElement& el, double J[3][3], double r, double s, double t);
    
    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjac0(const FESolidElement& el, double r, double s, double t, mat3d& J);
    
    //! calculate inverse jacobian matrix w.r.t. current frame
    double invjact(FESolidElement& el, double J[3][3], int n);
    
    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjact(FESolidElement& el, double J[3][3], double r, double s, double t);
    
    //! calculate inverse jacobian matrix w.r.t. current frame
    double invjact(FESolidElement& el, double J[3][3], int n, const double alpha);
    
    //! calculate inverse jacobian matrix w.r.t. previous time
    double invjactp(FESolidElement& el, double J[3][3], int n);
    
    //! calculate gradient of function at integration points
    vec3d gradient(FESolidElement& el, double* fn, int n);
    
    //! calculate gradient of function at integration points
    vec3d gradient(FESolidElement& el, vector<double>& fn, int n);
    
    //! calculate spatial gradient of vector function at integration points
    mat3d gradient(FESolidElement& el, vec3d* fn, int n);
    
    //! calculate spatial gradient of vector function at integration points
    //! at previous time
    mat3d gradientp(FESolidElement& el, vec3d* fn, int n);
    
    //! calculate material gradient of scalar function at integration points
    vec3d Gradient(FESolidElement& el, double* fn, int n);
    
    //! calculate material gradient of vector function at integration points
    mat3d Gradient(FESolidElement& el, vec3d* fn, int n);
    
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
    
    //! calculates contravariant basis vectors in reference configuration at an integration point
    void ContraBaseVectors0(FESolidElement& el, int j, vec3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FESolidElement& el, int j, vec3d g[3]);
    
    //! calculates parametric derivatives of covariant basis vectors at an integration point
    void CoBaseVectorDerivatives(FESolidElement& el, int j, vec3d dg[3][3]);
    
    //! calculates parametric derivatives of contravariant basis vectors at an integration point
    void ContraBaseVectorDerivatives(FESolidElement& el, int j, vec3d dg[3][3]);
    
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

	//! calculate the volume of an element
	double Volume(FESolidElement& el);

public:
	//! get the current nodal coordinates
	void GetCurrentNodalCoordinates(const FESolidElement& el, vec3d* rt);
	void GetCurrentNodalCoordinates(const FESolidElement& el, vec3d* rt, double alpha);

	//! get the reference nodal coordinates
	void GetReferenceNodalCoordinates(const FESolidElement& el, vec3d* r0);

	//! get the nodal coordinates at previous state
	void GetPreviousNodalCoordinates(const FESolidElement& el, vec3d* rp);

protected:
    vector<FESolidElement>	m_Elem;		//!< array of elements
    int     m_dofx;
    int     m_dofy;
    int     m_dofz;
    int     m_dofsx;
    int     m_dofsy;
    int     m_dofsz;
    int     m_dofsxp;
    int     m_dofsyp;
    int     m_dofszp;
};
