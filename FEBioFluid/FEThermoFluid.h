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
#include "FEElasticFluid.h"
#include "FEFluidThermalConductivity.h"
#include "FEFluidMaterialPoint.h"
#include "FEThermoFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
//! Base class for fluid materials.

class FEBIOFLUID_API FEThermoFluid : public FEFluidMaterial
{
public:
    FEThermoFluid(FEModel* pfem);
    
    // returns a pointer to a new material point object
    FEMaterialPointData* CreateMaterialPointData() override;
    
    //! Serialization
    void Serialize(DumpStream& ar) override;

public:
    //! calculate stress at material point
    mat3ds Stress(FEMaterialPoint& pt) override;
    
    //! tangent of stress with respect to strain J
    mat3ds Tangent_Strain(FEMaterialPoint& mp) override;
    
    //! tangent of stress with respect to temperature T
    mat3ds Tangent_Temperature(FEMaterialPoint& mp);
    
    //! Elastic pressure
    double Pressure(FEMaterialPoint& mp) override { return m_pElastic->Pressure(mp); }

    //! tangent of elastic pressure with respect to strain J
    double Tangent_Pressure_Strain(FEMaterialPoint& mp) override { return m_pElastic->Tangent_Strain(mp); }
    
    //! 2nd tangent of elastic pressure with respect to strain J
    double Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp) override { return m_pElastic->Tangent_Strain_Strain(mp); }
    
    //! tangent of elastic pressure with respect to temperature T
    double Tangent_Pressure_Temperature(FEMaterialPoint& mp) { return m_pElastic->Tangent_Temperature(mp); }
    
    //! 2nd tangent of elastic pressure with respect to temperature T
    double Tangent_Pressure_Temperature_Temperature(FEMaterialPoint& mp) { return m_pElastic->Tangent_Temperature_Temperature(mp); }
    
    //! tangent of elastic pressure with respect to strain J and temperature T
    double Tangent_Pressure_Strain_Temperature(FEMaterialPoint& mp) { return m_pElastic->Tangent_Strain_Temperature(mp); }
    
public:
    //! bulk modulus
    double BulkModulus(FEMaterialPoint& mp) override;

    //! heat flux
    vec3d HeatFlux(FEMaterialPoint& mp);
    
    //! strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp) override;
    
    //! invert pressure-dilatation relation
    bool Dilatation(const double T, const double p, double& e) override { return GetElastic()->Dilatation(T,p,e); }
    
    //! fluid pressure from state variables
    double Pressure(const double ef, const double T) override { return GetElastic()->Pressure(ef, T); };

    //! evaluate temperature
    double Temperature(FEMaterialPoint& mp) override;

public:
    //! return elastic part
    FEElasticFluid* GetElastic() { return m_pElastic; }
    
    //! return thermal conductivity part
    FEFluidThermalConductivity* GetConduct() { return m_pConduct; }
    
private: // material properties
    FEElasticFluid*             m_pElastic;     //!< pointer to elastic part of fluid material
    FEFluidThermalConductivity* m_pConduct;     //!< pointer to fluid thermal conductivity material
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
