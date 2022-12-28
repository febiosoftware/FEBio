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
#include "FEElasticFiberMaterialUC.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! D. E. Kiousis, T. C. Gasser and G. A. Holzapfel
//! Smooth contact strategies with emphasis on the modeling of balloon angioplasty with stenting
//! Int. J. Numer. Meth. Engng 2008; 75:826–855, equation (21)

class FEUncoupledFiberKiousis : public FEElasticFiberMaterialUC
{
public:
    FEUncoupledFiberKiousis(FEModel* pfem);
    
    //! Cauchy stress
    virtual mat3ds DevFiberStress(FEMaterialPoint& mp, const vec3d& a0) override;
    
    // Spatial tangent
    virtual tens4ds DevFiberTangent(FEMaterialPoint& mp, const vec3d& a0) override;
    
    //! Strain energy density
    virtual double DevFiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;
    
protected:
    FEParamDouble   m_d1;   // material coefficient d1 (proportional to initial fiber modulus)
    FEParamDouble   m_d2;   // material coefficient d2 (square of stretch when fiber engages)
    FEParamDouble   m_n;    // material coefficient n (power exponent)
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
};
