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



#include "FEFluidMaterial.h"
#include "FEFluidMaterialPoint.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEFluidMaterial, FEMaterial)

    // material parameters
    ADD_PARAMETER(m_rhor, FE_RANGE_GREATER_OR_EQUAL(0.0), "density")->setUnits(UNIT_DENSITY);

    // material properties
    ADD_PROPERTY(m_pViscous, "viscous");
//  EXPERIMENTAL
//    ADD_PROPERTY(m_pViscpol, "polar", FEProperty::Optional);

END_FECORE_CLASS();

//============================================================================
// FEFluid
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluid constructor

FEFluidMaterial::FEFluidMaterial(FEModel* pfem) : FEMaterial(pfem)
{
    m_rhor = 0;
    m_pViscous = nullptr;
    m_pViscpol = nullptr;
}

//-----------------------------------------------------------------------------
//! calculate current fluid density
double FEFluidMaterial::Density(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    return m_rhor/(vt.m_ef+1);
}

//-----------------------------------------------------------------------------
//! calculate current fluid kinematic viscosity
double FEFluidMaterial::KinematicViscosity(FEMaterialPoint& mp)
{
    return m_pViscous->ShearViscosity(mp)/Density(mp);
}

//-----------------------------------------------------------------------------
//! calculate current acoustic speed
double FEFluidMaterial::AcousticSpeed(FEMaterialPoint& mp)
{
    double c = sqrt(BulkModulus(mp)/Density(mp));
    
    return c;
}

//-----------------------------------------------------------------------------
//! calculate kinetic energy density (per reference volume)
double FEFluidMaterial::KineticEnergyDensity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double ked = m_rhor*(fp.m_vft*fp.m_vft)/2;
    return ked;
}

//-----------------------------------------------------------------------------
//! calculate strain + kinetic energy density (per reference volume)
double FEFluidMaterial::EnergyDensity(FEMaterialPoint& mp)
{
    return StrainEnergyDensity(mp) + KineticEnergyDensity(mp);
}

