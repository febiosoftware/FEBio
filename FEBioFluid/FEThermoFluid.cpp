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



#include "FEThermoFluid.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEThermoFluid, FEFluidMaterial)

    // material properties
    ADD_PROPERTY(m_pElastic, "elastic");
    ADD_PROPERTY(m_pConduct, "conduct");

END_FECORE_CLASS();

//============================================================================
// FEThermoFluid
//============================================================================

//-----------------------------------------------------------------------------
//! FEThermoFluid constructor

FEThermoFluid::FEThermoFluid(FEModel* pfem) : FEFluidMaterial(pfem)
{
    m_pElastic = 0;
    m_pConduct = 0;
}

//-----------------------------------------------------------------------------
void FEThermoFluid::Serialize(DumpStream& ar)
{
    FEFluidMaterial::Serialize(ar);
    if (ar.IsShallow()) return;
}

//-----------------------------------------------------------------------------
//! returns a pointer to a new material point object
FEMaterialPointData* FEThermoFluid::CreateMaterialPointData()
{
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    return new FEThermoFluidMaterialPoint(fp);
}

//-----------------------------------------------------------------------------
//! evaluate temperature
double FEThermoFluid::Temperature(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tp = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    return tp.m_T;
}

//-----------------------------------------------------------------------------
//! bulk modulus
double FEThermoFluid::BulkModulus(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    return -(vt.m_ef+1)*Tangent_Pressure_Strain(mp);
}

//-----------------------------------------------------------------------------
//! heat flux
vec3d FEThermoFluid::HeatFlux(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tp = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double k = m_pConduct->ThermalConductivity(mp);
    vec3d q = -tp.m_gradT*k;
    return q;
}

//-----------------------------------------------------------------------------
//! The stress of a fluid material is the sum of the fluid pressure
//! and the viscous stress.

mat3ds FEThermoFluid::Stress(FEMaterialPoint& mp)
{
    // calculate solid material stress
    mat3ds s = GetViscous()->Stress(mp);
    
    double p = m_pElastic->Pressure(mp);
    
    // add fluid pressure
    s.xx() -= p;
    s.yy() -= p;
    s.zz() -= p;
    
    return s;
}

//-----------------------------------------------------------------------------
//! The tangent of stress with respect to strain J of a fluid material is the
//! sum of the tangent of the fluid pressure and that of the viscous stress.

mat3ds FEThermoFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    // get tangent of viscous stress
    mat3ds sJ = GetViscous()->Tangent_Strain(mp);
    
    // add tangent of fluid pressure
    double dp = m_pElastic->Tangent_Strain(mp);
    sJ.xx() -= dp;
    sJ.yy() -= dp;
    sJ.zz() -= dp;
    
    return sJ;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density (per reference volume)
double FEThermoFluid::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = m_rhor*m_pElastic->SpecificStrainEnergy(mp);
    return sed;
}
