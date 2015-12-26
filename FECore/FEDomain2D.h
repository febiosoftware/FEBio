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
class FEDomain2D : public FEDomain
{
public:
    //! constructor
    FEDomain2D(FEMesh* pm) : FEDomain(FE_DOMAIN_2D, pm) {}
    
    //! create storage for elements
    void create(int nsize) { m_Elem.resize(nsize); }
    
    //! return nr of elements
    int Elements() { return m_Elem.size(); }
    
    //! element access
    FEElement2D& Element(int n) { return m_Elem[n]; }
    FEElement& ElementRef(int n) { return m_Elem[n]; }
    
    int GetElementType() { return m_Elem[0].Type(); }
    
    //! Initialize elements
    void InitElements();
    
    //! Reset element data
    void Reset();
    
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
    
public:
    //! shallow copy
    void ShallowCopy(DumpStream& dmp, bool bsave);
    
    //! Serialize domain data to archive
    void Serialize(DumpFile& ar);
    
protected:
    vector<FEElement2D>	m_Elem;	//!< array of elements
};
