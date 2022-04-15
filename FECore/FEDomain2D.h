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

//-----------------------------------------------------------------------------
//! Abstract base class for shell elements
class FECORE_API FEDomain2D : public FEDomain
{
    FECORE_SUPER_CLASS(FEDOMAIN2D_ID)
    FECORE_BASE_CLASS(FEDomain2D)

public:
    //! constructor
    FEDomain2D(FEModel* fem) : FEDomain(FE_DOMAIN_2D, fem) {}
    
    //! create storage for elements
	bool Create(int nsize, FE_Element_Spec espec) override;

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
