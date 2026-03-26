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
#include "FENewtonianThermoFluid.h"
#include "FEThermoFluid.h"
#include "FEThermalPropConst.h"
#include <FEBioFluid/FEFluid.h>
#include <FEBioFluid/FEBiphasicFSI.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>

// define the material parameters
BEGIN_FECORE_CLASS(FENewtonianThermoFluid, FEThermoViscousFluid)
    ADD_PARAMETER(m_kappa, FE_RANGE_GREATER_OR_EQUAL(0.0), "kappa")->setUnits(UNIT_VISCOSITY)->setLongName("referential bulk viscosity");
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"   )->setUnits(UNIT_VISCOSITY)->setLongName("referential shear viscosity");
// Optionally add strain-dependent normalized viscosity relations
    ADD_PROPERTY(m_kappahat, "kappahat", FEProperty::Optional)->SetLongName("normalized bulk viscosity");
    ADD_PROPERTY(m_muhat, "muhat", FEProperty::Optional)->SetLongName("normalized bulk viscosity");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FENewtonianThermoFluid::FENewtonianThermoFluid(FEModel* pfem) : FEThermoViscousFluid(pfem)
{
    m_kappahat = nullptr;
    m_muhat = nullptr;
}

//-----------------------------------------------------------------------------
//! initialization
bool FENewtonianThermoFluid::Init()
{
    FEModel* pfem = GetFEModel();
    if (m_kappahat == nullptr) {
        m_kappahat = new FEThermalPropConst(pfem);
    }
    m_kappahat->Init();
    if (m_muhat == nullptr) {
        m_muhat = new FEThermalPropConst(pfem);
    }
    m_muhat->Init();
    m_Tr = GetFEModel()->GetGlobalConstant("T");
    if (m_Tr <= 0) { feLogError("A positive absolute temperature T must be defined in Globals section");     return false; }

    return FEViscousFluid::Init();
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FENewtonianThermoFluid::Stress(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& vt = *pt.ExtractData<FEFluidMaterialPoint>();
    
    mat3ds D = vt.RateOfDeformation();
    
    double mu = ShearViscosity(pt);
    double kappa = BulkViscosity(pt);
    
    mat3ds s = mat3dd(1.0)*(D.tr()*(kappa - 2.*mu/3.)) + D*(2*mu);
        
    return s;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to strain J
mat3ds FENewtonianThermoFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    
    mat3ds D = vt.RateOfDeformation();
    
    double dmudJ = m_mu*m_muhat->Tangent_NormalizedProperty_Strain(mp);
    double dkappadJ = m_kappa*m_kappahat->Tangent_NormalizedProperty_Strain(mp);
    
    mat3ds dsdJ = mat3dd(1.0)*(D.tr()*(dkappadJ - 2.*dmudJ/3.)) + D*(2*dmudJ);
        
    return dsdJ;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FENewtonianThermoFluid::Tangent_RateOfDeformation(FEMaterialPoint& mp)
{
    mat3dd I(1.0);
    double mu = ShearViscosity(mp);
    double kappa = BulkViscosity(mp);
    
    tens4ds c = dyad1s(I)*(kappa - 2.*mu/3.) + dyad4s(I)*(2*mu);
    
    return c;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to temperature T
mat3ds FENewtonianThermoFluid::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    
    mat3ds D = vt.RateOfDeformation();
    
    double dkappa = 0;
    double dmu = 0;
    
    mat3ds ds = mat3dd(1.0)*(D.tr()*(dkappa - 2.*dmu/3.)) + D*(2*dmu);
        
    return ds;
}

//-----------------------------------------------------------------------------
//! calculate viscosities
double FENewtonianThermoFluid::ShearViscosity(FEMaterialPoint& mp) {
    return m_mu*m_muhat->NormalizedProperty(mp);
}

double FENewtonianThermoFluid::TangentShearViscosityTemperature(FEMaterialPoint& mp) {
    return m_mu*m_muhat->Tangent_NormalizedProperty_Temperature(mp);
}

double FENewtonianThermoFluid::TangentShearViscosityStrain(FEMaterialPoint& mp) {
    return m_mu*m_muhat->Tangent_NormalizedProperty_Strain(mp);
}

double FENewtonianThermoFluid::BulkViscosity(FEMaterialPoint& mp) {
    return m_kappa*m_muhat->NormalizedProperty(mp);
}

double FENewtonianThermoFluid::TangentBulkViscosityTemperature(FEMaterialPoint& mp) {
    return m_kappa*m_muhat->Tangent_NormalizedProperty_Temperature(mp);
}

double FENewtonianThermoFluid::TangentBulkViscosityStrain(FEMaterialPoint& mp) {
    return m_kappa*m_muhat->Tangent_NormalizedProperty_Strain(mp);
}
