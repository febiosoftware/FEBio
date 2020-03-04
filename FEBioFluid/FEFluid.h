/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include <FECore/FEModel.h>
#include <FEBioMech/FEBodyForce.h>
#include "FEFluidMaterial.h"
#include "FEFluidMaterialPoint.h"
#include "FEViscousFluid.h"

//-----------------------------------------------------------------------------
//! Base class for fluid materials.

class FEBIOFLUID_API FEFluid : public FEFluidMaterial
{
public:
	FEFluid(FEModel* pfem);
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() override;
	
public:
    //! initialization
    bool Init() override {
        m_Tr = GetFEModel()->GetGlobalConstant("T");
        return true;
    }
    
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
	
    //! tangent of stress with respect to strain J
    mat3ds Tangent_Strain(FEMaterialPoint& mp) override;
    
    //! elastic pressure
    double Pressure(FEMaterialPoint& mp) override;
    virtual double Pressure(const double e) override;

    //! tangent of elastic pressure with respect to strain J
    double Tangent_Pressure_Strain(FEMaterialPoint& mp) override { return -m_k; }
    
    //! 2nd tangent of elastic pressure with respect to strain J
    double Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp) override { return 0; }
    
    //! bulk modulus
    double BulkModulus(FEMaterialPoint& mp) override;
    
    //! strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp) override;
    
    //! invert pressure-dilatation relation
    double Dilatation(const double p) override;
    
    //! evaluate temperature
    double Temperature(FEMaterialPoint& mp) override { return m_Tr; }
    
public:
    double      m_k;        //!< bulk modulus at J=1
    double      m_Tr;       //!< ambient temperature
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
