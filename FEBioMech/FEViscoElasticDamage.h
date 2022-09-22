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
#include "FEViscoElasticMaterial.h"
#include "FEDamageMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a large deformation visco-elastic material
//
class FEViscoElasticDamage : public FEElasticMaterial
{
public:
    // NOTE: make sure that this parameter is the
    //       same as the MAX_TERMS in the FEViscoElasticMaterialPoint class
    enum { MAX_TERMS = FEViscoElasticMaterialPoint::MAX_TERMS };
    
public:
    //! default constructor
    FEViscoElasticDamage(FEModel* pfem);
    
public:
    //! initialization
    bool Init() override;
    
    //! stress function
    mat3ds Stress(FEMaterialPoint& pt) override;
    
    //! tangent function
    tens4ds Tangent(FEMaterialPoint& pt) override;
    
    //! strain energy density
    double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! calculate exponent of right-stretch tensor in series spring
    bool SeriesStretchExponent(FEMaterialPoint& pt);
    
    // returns a pointer to a new material point object
    FEMaterialPointData* CreateMaterialPointData() override;
    
public:
    // material parameters
    double    m_g0;            //!< intitial visco-elastic coefficient
    double    m_g[MAX_TERMS];    //!< visco-elastic coefficients
    double    m_t[MAX_TERMS];    //!< relaxation times
    
private:
    FEDamageMaterial*    m_pDmg;    //!< pointer to elastic damage material

public:
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
