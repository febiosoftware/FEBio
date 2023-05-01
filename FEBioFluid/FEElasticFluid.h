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
#include <FECore/FEMaterial.h>
#include <FECore/tens4d.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! Base class for the viscous part of the fluid response.
//! These materials provide the viscous stress and its tangents.
//!
class FEBIOFLUID_API FEElasticFluid : public FEMaterialProperty
{
public:
    FEElasticFluid(FEModel* pfem) : FEMaterialProperty(pfem) {}
    virtual ~FEElasticFluid() {}
    
    //! gage pressure
    virtual double Pressure(FEMaterialPoint& pt) = 0;
    
    //! tangent of pressure with respect to strain J
    virtual double Tangent_Strain(FEMaterialPoint& mp) = 0;
    
    //! 2nd tangent of pressure with respect to strain J
    virtual double Tangent_Strain_Strain(FEMaterialPoint& mp) = 0;
    
    //! tangent of pressure with respect to temperature T
    virtual double Tangent_Temperature(FEMaterialPoint& mp) = 0;
    
    //! 2nd tangent of pressure with respect to temperature T
    virtual double Tangent_Temperature_Temperature(FEMaterialPoint& mp) = 0;
    
    //! tangent of pressure with respect to strain J and temperature T
    virtual double Tangent_Strain_Temperature(FEMaterialPoint& mp) = 0;
    
    //! specific free energy
    virtual double SpecificFreeEnergy(FEMaterialPoint& mp) = 0;
    
    //! specific entropy
    virtual double SpecificEntropy(FEMaterialPoint& mp) = 0;
    
    //! specific strain energy
    virtual double SpecificStrainEnergy(FEMaterialPoint& mp) = 0;
    
    //! isochoric specific heat capacity
    virtual double IsochoricSpecificHeatCapacity(FEMaterialPoint& mp) = 0;
            
    //! tangent of isochoric specific heat capacity with respect to strain J
    virtual double Tangent_cv_Strain(FEMaterialPoint& mp) = 0;
            
    //! tangent of isochoric specific heat capacity with respect to temperature T
    virtual double Tangent_cv_Temperature(FEMaterialPoint& mp) = 0;

    //! isobaric specific heat capacity
    virtual double IsobaricSpecificHeatCapacity(FEMaterialPoint& mp) = 0;
            
    //! calculate dilatation for given (effective) pressure and temperature
    virtual bool Dilatation(const double T, const double p, double& e) = 0;
    
    //! calculate fluid pressure and its derivatives from state variables
    double Pressure(const double ef, const double T);
    double Tangent_Strain(const double ef, const double T);
    double Tangent_Temperature(const double ef, const double T);
    double Tangent_Strain_Strain(const double ef, const double T);
    double Tangent_Strain_Temperature(const double ef, const double T);
    double Tangent_Temperature_Temperature(const double ef, const double T);

public:
    //! specific internal energy
    double SpecificInternalEnergy(FEMaterialPoint& mp);
    
    //! specific gage enthalpy
    double SpecificGageEnthalpy(FEMaterialPoint& mp);
    
    //! specific free enthalpy
    double SpecificFreeEnthalpy(FEMaterialPoint& mp);

    FECORE_BASE_CLASS(FEElasticFluid)
};
