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
#include "FEElasticFluid.h"
#include <FECore/FEFunction1D.h>

//-----------------------------------------------------------------------------
//! Linear elastic fluid
//!
class FEBIOFLUID_API FELogNonlinearElasticFluid : public FEElasticFluid
{
public:
    FELogNonlinearElasticFluid(FEModel* pfem);
    ~FELogNonlinearElasticFluid() {}
    
    //! initialization
    bool Init() override;
    
    // serialization
    void Serialize(DumpStream& ar) override;

    //! gage pressure
    double Pressure(FEMaterialPoint& pt) override;
    
    //! tangent of pressure with respect to strain J
    double Tangent_Strain(FEMaterialPoint& mp) override;
    
    //! 2nd tangent of pressure with respect to strain J
    double Tangent_Strain_Strain(FEMaterialPoint& mp) override;
    
    //! tangent of pressure with respect to temperature T
    double Tangent_Temperature(FEMaterialPoint& mp) override;
    
    //! 2nd tangent of pressure with respect to temperature T
    double Tangent_Temperature_Temperature(FEMaterialPoint& mp) override;
    
    //! tangent of pressure with respect to strain J and temperature T
    double Tangent_Strain_Temperature(FEMaterialPoint& mp) override;
    
    //! specific free energy
    double SpecificFreeEnergy(FEMaterialPoint& mp) override;
    
    //! specific entropy
    double SpecificEntropy(FEMaterialPoint& mp) override;
    
    //! specific strain energy
    double SpecificStrainEnergy(FEMaterialPoint& mp) override;
    
    //! isochoric specific heat capacity
    double IsochoricSpecificHeatCapacity(FEMaterialPoint& mp) override;
    
    //! tangent of isochoric specific heat capacity with respect to strain J
    double Tangent_cv_Strain(FEMaterialPoint& mp) override;
    
    //! tangent of isochoric specific heat capacity with respect to temperature T
    double Tangent_cv_Temperature(FEMaterialPoint& mp) override;
    
    //! isobaric specific heat capacity
    double IsobaricSpecificHeatCapacity(FEMaterialPoint& mp) override;
    
    //! dilatation from temperature and pressure
    bool Dilatation(const double T, const double p, double& e) override;
    
public:
    double      m_k;        //!< bulk modulus
    double      m_rhor;     //!< true density
};
