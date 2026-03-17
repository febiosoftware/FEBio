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
#include "FEThermoFluidMaterial.h"
#include "FEThermoFluidMaterialPoint.h"
#include "FEThermoViscousFluid.h"
#include <FEBioFluid/FEFluidMaterialPoint.h>
#include "FEThermoElasticFluid.h"

//-----------------------------------------------------------------------------
//! Base class for thermo-fluid materials.

class FEBIOTHERMOFLUID_API FEThermoFluid : public FEThermoFluidMaterial
{
public:
    FEThermoFluid(FEModel* pfem);
	
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
    
    //! tangent of elastic pressure with respect to temperature T
    double Tangent_Pressure_Temperature(FEMaterialPoint& mp) { return m_pElastic->Tangent_Temperature(mp); }
    
    //! 2nd tangent of elastic pressure with respect to temperature T
    double Tangent_Pressure_Temperature_Temperature(FEMaterialPoint& mp) { return m_pElastic->Tangent_Temperature_Temperature(mp); }
    
    //! tangent of elastic pressure with respect to strain J and temperature T
    double Tangent_Pressure_Strain_Temperature(FEMaterialPoint& mp) { return m_pElastic->Tangent_Strain_Temperature(mp); }
    
public:
    //! bulk modulus
    double BulkModulus(FEMaterialPoint& mp) override;
    
    //! strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp) override;
    
    //! evaluate dilatation from effective pressure
    bool Dilatation(const double T, const double p, double& e) override;
    
    //! return elastic fluid
    FEThermoElasticFluid* GetElastic() { return m_pElastic; }
    
private: // material properties
    FEThermoElasticFluid*           m_pElastic;     //!< pointer to elastic part of fluid material

public:
    double      m_Tr;       //!< ambient temperature
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
