/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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
#include <FECore/FEMaterial.h>
#include <FECore/tens4d.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! Base class for the viscous part of the fluid response.
//! These materials provide the viscous stress and its tangents.
//!
class FEBIOFLUID_API FEViscousPolarFluid : public FEMaterial
{
public:
    FEViscousPolarFluid(FEModel* pfem) : FEMaterial(pfem) {}
    virtual ~FEViscousPolarFluid() {}
    
    //! dual vector non-symmetric part of viscous stress in polar fluid
    virtual vec3d SkewStressDualVector(FEMaterialPoint& pt) = 0;
    
    //! non-symmetric part of viscous stress in polar fluid
    virtual mat3da SkewStress(FEMaterialPoint& pt) { return mat3da(SkewStressDualVector(pt)); }
    
    //! tangent of stress with respect to strain J
    virtual mat3da SkewTangent_Strain(FEMaterialPoint& mp) = 0;
    
    //! tangent of stress with respect to relative rate of rotation tensor H
    virtual mat3d SkewTangent_RateOfRotation(FEMaterialPoint& mp) = 0;
    
    //! tangent of stress with respect to temperature
    virtual mat3da SkewTangent_Temperature(FEMaterialPoint& mp) = 0;
    
    //! viscous couple stress in polar fluid
    virtual mat3d CoupleStress(FEMaterialPoint& pt) = 0;
    
    //! tangent of viscous couple stress with respect to strain J
    virtual mat3d CoupleTangent_Strain(FEMaterialPoint& mp) = 0;
    
    //! tangent of viscous couple stress with respect to rate of rotation tensor g
    virtual tens4d CoupleTangent_RateOfRotation(FEMaterialPoint& mp) = 0;
    
    //! tangent of viscous couple stress with respect to temperature
    virtual mat3d CoupleTangent_Temperature(FEMaterialPoint& mp) = 0;
    
    //! dynamic viscosity
    virtual double RelativeRotationalViscosity(FEMaterialPoint& mp) = 0;
    
};
