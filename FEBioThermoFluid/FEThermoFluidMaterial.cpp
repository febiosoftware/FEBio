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
#include "FEThermoFluidMaterial.h"
#include "FEThermoFluidMaterialPoint.h"
#include <FEBioFluid/FEFluidMaterialPoint.h>
#include <FECore/DumpStream.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEThermoFluidMaterial, FEMaterial)

// material parameters
    ADD_PARAMETER(m_rhor, FE_RANGE_GREATER_OR_EQUAL(0.0), "density")->setUnits(UNIT_DENSITY);


// material properties
    ADD_PROPERTY(m_pViscous, "viscous");
    ADD_PROPERTY(m_pHeatFlux, "heat_flux");

END_FECORE_CLASS();

//============================================================================
// FEThermoFluidMaterial
//============================================================================

//-----------------------------------------------------------------------------
//! FEThermoFluidMaterial constructor

FEThermoFluidMaterial::FEThermoFluidMaterial(FEModel* pfem) : FEMaterial(pfem)
{
    m_rhor = 0;
    m_pHeatFlux = nullptr;
    m_pViscous = nullptr;
}

//-----------------------------------------------------------------------------
//! evaluate temperature
double FEThermoFluidMaterial::Temperature(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tp = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    return tp.m_T;
}

//-----------------------------------------------------------------------------
//! bulk modulus
double FEThermoFluidMaterial::BulkModulus(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    return -(vt.m_ef+1)*Tangent_Pressure_Strain(mp);
}

//-----------------------------------------------------------------------------
//! heat flux
vec3d FEThermoFluidMaterial::HeatFlux(FEMaterialPoint& mp)
{
    return m_pHeatFlux->HeatFlux(mp);;
}

//-----------------------------------------------------------------------------
//! calculate current fluid density
double FEThermoFluidMaterial::Density(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    return m_rhor/(vt.m_ef+1);
}

//-----------------------------------------------------------------------------
//! calculate current fluid kinematic viscosity
double FEThermoFluidMaterial::KinematicViscosity(FEMaterialPoint& mp)
{
    return m_pViscous->ShearViscosity(mp)/Density(mp);
}

//-----------------------------------------------------------------------------
//! calculate current acoustic speed
double FEThermoFluidMaterial::AcousticSpeed(FEMaterialPoint& mp)
{
    double c = sqrt(BulkModulus(mp)/Density(mp));
    
    return c;
}

//-----------------------------------------------------------------------------
//! calculate kinetic energy density (per reference volume)
double FEThermoFluidMaterial::KineticEnergyDensity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double ked = m_rhor*(fp.m_vft*fp.m_vft)/2;
    return ked;
}

//-----------------------------------------------------------------------------
//! calculate strain + kinetic energy density (per reference volume)
double FEThermoFluidMaterial::EnergyDensity(FEMaterialPoint& mp)
{
    return StrainEnergyDensity(mp) + KineticEnergyDensity(mp);
}

