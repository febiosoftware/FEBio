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
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
//! A material class describing the active fiber contraction model of Estrada et al. in doi: 10.1115/1.4044030
class FEBIOMECH_API FEActiveContractionMaterial : public FEMaterialProperty
{
public:
    FEActiveContractionMaterial(FEModel* pfem) : FEMaterialProperty(pfem) {}
    virtual ~FEActiveContractionMaterial(){}
    
    //! calculate the active contractile stress
    virtual mat3ds ActiveStress(FEMaterialPoint& mp, const vec3d& a0) = 0;
    
    //! active contraction stiffness contribution
    virtual tens4ds ActiveStiffness(FEMaterialPoint& mp, const vec3d& a0) = 0;
    
    //! update force-velocity material point
    virtual void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp, const vec3d& a0) {}

    FECORE_BASE_CLASS(FEActiveContractionMaterial)
};
