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
#include "FEViscousPolarFluid.h"

//-----------------------------------------------------------------------------
// This class evaluates the viscous stress in a linear polar fluid

class FEBIOFLUID_API FEViscousPolarLinear : public FEViscousPolarFluid
{
public:
    //! constructor
    FEViscousPolarLinear(FEModel* pfem);
    
    //! dual vector non-symmetric part of viscous stress in polar fluid
    vec3d SkewStressDualVector(FEMaterialPoint& pt) override;
    
    //! non-symmetric part of viscous stress in polar fluid
    mat3da SkewStress(FEMaterialPoint& pt) override;
    
    //! tangent of stress with respect to strain J
    mat3da SkewTangent_Strain(FEMaterialPoint& mp) override;
    
    //! tangent of stress with respect to relative rate of rotation tensor H
    mat3d SkewTangent_RateOfRotation(FEMaterialPoint& mp) override;
    
    //! tangent of stress with respect to temperature
    mat3da SkewTangent_Temperature(FEMaterialPoint& mp) override { return mat3da(0,0,0); }
    
    //! viscous couple stress in polar fluid
    mat3d CoupleStress(FEMaterialPoint& pt) override;
    
    //! tangent of viscous couple stress with respect to strain J
    mat3d CoupleTangent_Strain(FEMaterialPoint& mp) override;
    
    //! tangent of viscous couple stress with respect to rate of rotation tensor g
    tens4d CoupleTangent_RateOfRotation(FEMaterialPoint& mp) override;
    
    //! tangent of viscous couple stress with respect to temperature
    mat3d CoupleTangent_Temperature(FEMaterialPoint& mp) override { return mat3d(0.); }
    
    //! dynamic viscosity
    double RelativeRotationalViscosity(FEMaterialPoint& mp) override;
    
public:
    double    m_tau;    //!< relative rotational viscosity
    double    m_alpha;  //!< couple stress viscosity
    double    m_beta;   //!< couple stress viscosity
    double    m_gamma;  //!< couple stress viscosity

    // declare parameter list
    DECLARE_FECORE_CLASS();
};
