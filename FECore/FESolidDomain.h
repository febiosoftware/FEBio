#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! abstract base class for 3D volumetric elements
class FESolidDomain : public FEDomain
{
public:
	//! constructor
	FESolidDomain(FEMesh* pm) : FEDomain(FE_DOMAIN_SOLID, pm) {}

	//! create storage for elements
	void create(int nsize) { m_Elem.resize(nsize); }

	//! return nr of elements
	int Elements() { return m_Elem.size(); }

	//! reset data (overridden from FEDomain)
	void Reset();

	//! copy data from another domain (overridden from FEDomain)
	void CopyFrom(FEDomain* pd);

	//! element access
	FESolidElement& Element(int n) { return m_Elem[n]; }
	FEElement& ElementRef(int n) { return m_Elem[n]; }

	int GetElementType() { return m_Elem[0].Type(); }

	//! find the element in which point y lies
	FESolidElement* FindElement(vec3d y, double r[3]);

	//! Calculate deformation gradient at integration point n
	double defgrad(FESolidElement& el, mat3d& F, int n);

	//! Calculate deformation gradient at integration point n
	double defgrad(FESolidElement& el, mat3d& F, double r, double s, double t);

	//! calculate inverse jacobian matrix w.r.t. reference frame
	double invjac0(FESolidElement& el, double J[3][3], int n);

	//! calculate inverse jacobian matrix w.r.t. reference frame
	double invjac0(const FESolidElement& el, double J[3][3], double r, double s, double t);

	//! calculate inverse jacobian matrix w.r.t. reference frame
	double invjac0(const FESolidElement& el, double r, double s, double t, mat3d& J);

	//! calculate inverse jacobian matrix w.r.t. current frame
	double invjact(FESolidElement& el, double J[3][3], int n);

	//! calculate inverse jacobian matrix w.r.t. reference frame
	double invjact(FESolidElement& el, double J[3][3], double r, double s, double t);

	//! calculate gradient of function at integration points
	vec3d gradient(FESolidElement& el, double* fn, int n);

	//! calculate gradient of function at integration points
	vec3d gradient(FESolidElement& el, vector<double>& fn, int n);

    //! calculate spatial gradient of vector function at integration points
    mat3d gradient(FESolidElement& el, vec3d* fn, int n);
    
    //! calculate material gradient of vector function at integration points
    mat3d Gradient(FESolidElement& el, vec3d* fn, int n);
    
	//! calculate jacobian in reference frame
	double detJ0(FESolidElement& el, int n);

	//! calculate jacobian in current frame
	double detJt(FESolidElement& el, int n);
    
    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FESolidElement& el, int j, vec3d g[3]);

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
    
public:
	//! serialize data to archive
	void Serialize(DumpStream& ar);

protected:
	vector<FESolidElement>	m_Elem;		//!< array of elements
};
