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
#include <FEBioMech/FEBodyForce.h>
#include "FEFluidMaterial.h"
#include "FEFluidMaterialPoint.h"
#include "FEViscousFluid.h"
#include <FEBioMix/FEBiphasic.h>
#include "FEElasticFluid.h"

//-----------------------------------------------------------------------------
//! Base class for fluid materials.

//! NOTE: This inherits from FEBiphasicInterface in order to override the GetActualFluidPressure, 
//!       which is used in FEReactionRateExpSED and FEReactionRateHuiskes. 
//!       Note sure yet if there is a better alternative.

class FEBIOFLUID_API FEFluid : public FEFluidMaterial, public FEBiphasicInterface
{
public:
	FEFluid(FEModel* pfem);
	
	// returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;
	
public:
    //! initialization
    bool Init() override;
    
    //! Serialization
    void Serialize(DumpStream& ar) override;

	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
	
    //! tangent of stress with respect to strain J
    mat3ds Tangent_Strain(FEMaterialPoint& mp) override;
    
    //! elastic pressure
    double Pressure(FEMaterialPoint& mp) override;
    double Pressure(const double e, const double T = 0) override;

    //! tangent of elastic pressure with respect to strain J
    double Tangent_Pressure_Strain(FEMaterialPoint& mp) override;
    
    //! 2nd tangent of elastic pressure with respect to strain J
    double Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp) override;
    
    //! bulk modulus
    double BulkModulus(FEMaterialPoint& mp) override;
    
    //! strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp) override;
    
    //! evaluate temperature
    double Temperature(FEMaterialPoint& mp) override { return m_Tr; }

    //! evaluate dilatation from effective pressure
    bool Dilatation(const double T, const double p, double& e) override;
    
    //! return elastic fluid
    FEElasticFluid* GetElastic() { return m_pElastic; }
    
private: // material properties
    FEElasticFluid*             m_pElastic;     //!< pointer to elastic part of fluid material
    

public: // from FEBiphasicInterface
    double GetActualFluidPressure(const FEMaterialPoint& mp) override {
        const FEFluidMaterialPoint* pt = (mp.ExtractData<FEFluidMaterialPoint>());
        return pt->m_pf;
    }
    
public:
    double      m_k;        //!< bulk modulus at J=1
    double      m_Tr;       //!< ambient temperature
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
