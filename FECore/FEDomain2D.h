//
//  FEDomain2D.hpp
//  FECore
//
//  Created by Gerard Ateshian on 12/17/15.
//  Copyright © 2015 febio.org. All rights reserved.
//

#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Abstract base class for shell elements
class FECORE_API FEDomain2D : public FEDomain
{
public:
    //! constructor
    FEDomain2D(FEModel* fem) : FEDomain(FE_DOMAIN_2D, fem) {}
    
    //! create storage for elements
    void Create(int nelems, int elemType) override;
    
    //! return nr of elements
    int Elements() const override { return (int)m_Elem.size(); }
    
    //! element access
    FEElement2D& Element(int n) { return m_Elem[n]; }
    FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

    int GetElementType() { return m_Elem[0].Type(); }
    
    //! Initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
    //! Reset element data
    void Reset() override;
    
    // inverse jacobian with respect to reference frame
    double invjac0(FEElement2D& el, double J[2][2], int n);
    
    // inverse jacobian with respect to current frame
    double invjact(FEElement2D& el, double J[2][2], int n);
    
    //! calculate in-plane gradient of function at integration points
    vec2d gradient(FEElement2D& el, double* fn, int n);
    
    //! calculate in-plane gradient of function at integration points
    vec2d gradient(FEElement2D& el, vector<double>& fn, int n);

    //! calculate in-plane gradient of vector function at integration points
    mat2d gradient(FEElement2D& el, vec2d* fn, int n);
    
    //! calculate in-plane gradient of vector function at integration points
    mat3d gradient(FEElement2D& el, vec3d* fn, int n);
    
    // jacobian with respect to reference frame
    double detJ0(FEElement2D& el, int n);
    
    //! calculate jacobian in current frame
    double detJt(FEElement2D& el, int n);
    
    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FEElement2D& el, int j, vec2d g[2]);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FEElement2D& el, int j, vec2d g[2]);
    
    //! calculates parametric derivatives of covariant basis vectors at an integration point
    void CoBaseVectorDerivatives(FEElement2D& el, int j, vec2d dg[2][2]);
    
    //! calculates parametric derivatives of contravariant basis vectors at an integration point
    void ContraBaseVectorDerivatives(FEElement2D& el, int j, vec2d dg[2][2]);
    
    //! calculate the laplacian of a vector function at an integration point
    vec2d lapvec(FEElement2D& el, vec2d* fn, int n);
    
    //! calculate the gradient of the divergence of a vector function at an integration point
    vec2d gradivec(FEElement2D& el, vec2d* fn, int n);
    
protected:
    vector<FEElement2D>	m_Elem;	//!< array of elements
};
