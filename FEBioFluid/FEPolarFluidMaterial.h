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
#include "FEFluidMaterial.h"
#include <FEBioMech/FEBodyForce.h>
#include "FEFluidMaterialPoint.h"
#include "FEViscousFluid.h"
#include "FEViscousPolarFluid.h"

//-----------------------------------------------------------------------------
//! Base class for fluid materials.

class FEBIOFLUID_API FEPolarFluidMaterial : public FEFluidMaterial
{
public:
    FEPolarFluidMaterial(FEModel* pfem);
    virtual ~FEPolarFluidMaterial() {}
    
public:
    //! calculate stress at material point
    virtual mat3ds Stress(FEMaterialPoint& pt) override = 0;
    
    //! tangent of stress with respect to strain J
    virtual mat3ds Tangent_Strain(FEMaterialPoint& mp) override = 0;
    
    //! elastic fluid pressure
    virtual double Pressure(FEMaterialPoint& mp) override = 0;
    
    //! tangent of elastic pressure with respect to strain J
    virtual double Tangent_Pressure_Strain(FEMaterialPoint& mp) override = 0;
    
    //! 2nd tangent of elastic pressure with respect to strain J
    virtual double Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp) override = 0;
    
    //! bulk modulus
    virtual double BulkModulus(FEMaterialPoint& mp) override = 0;
    
    //! strain energy density
    virtual double StrainEnergyDensity(FEMaterialPoint& mp) override = 0;
    
    //! invert effective pressure-dilatation relation
    virtual bool Dilatation(const double T, const double p, double& e) override = 0;
    
    //! evaluate temperature
    virtual double Temperature(FEMaterialPoint& mp) override = 0;
    
    //! return viscous part
    FEViscousPolarFluid* GetViscousPolar() { return m_pViscpol; }
    
    //! fluid pressure from state variables
    virtual double Pressure(const double ef, const double T) override = 0;
    
private: // material properties
    FEViscousPolarFluid*    m_pViscpol; //!< pointer to viscous polar part of fluid material

    // declare parameter list
    DECLARE_FECORE_CLASS();
};
