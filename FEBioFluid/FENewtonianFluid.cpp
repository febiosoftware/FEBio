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



#include "stdafx.h"
#include "FENewtonianFluid.h"
#include "FEFluid.h"
#include "FEBiphasicFSI.h"
#include "FEThermoFluid.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FENewtonianFluid, FEViscousFluid)
    ADD_PARAMETER(m_kappa, FE_RANGE_GREATER_OR_EQUAL(0.0), "kappa")->setUnits(UNIT_VISCOSITY)->setLongName("bulk viscosity");
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"   )->setUnits(UNIT_VISCOSITY)->setLongName("shear viscosity");
// properties
    ADD_PROPERTY(m_kappahat, "khat" ,FEProperty::Optional)->SetLongName("normalized bulk viscosity");
    ADD_PROPERTY(m_muhat   , "muhat",FEProperty::Optional)->SetLongName("normalized shear viscosity");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FENewtonianFluid::FENewtonianFluid(FEModel* pfem) : FEViscousFluid(pfem)
{
    m_kappa = 0;
    m_mu = 0;
    m_Tr = 0;
    m_kappahat = nullptr;
    m_muhat = nullptr;
}

//-----------------------------------------------------------------------------
//! initialization
bool FENewtonianFluid::Init()
{
    if (m_kappahat || m_muhat) {
        m_Tr = GetGlobalConstant("T");
        
        if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    }
    
    if (m_kappahat) m_kappahat->Init();
    if (m_muhat) m_muhat->Init();
    
    return FEViscousFluid::Init();
}

//-----------------------------------------------------------------------------
void FENewtonianFluid::Serialize(DumpStream& ar)
{
    FEViscousFluid::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    ar & m_Tr;
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FENewtonianFluid::Stress(FEMaterialPoint& pt)
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
mat3ds FENewtonianFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return mat3ds(0,0,0,0,0,0);
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FENewtonianFluid::Tangent_RateOfDeformation(FEMaterialPoint& mp)
{
    mat3dd I(1.0);
    double mu = ShearViscosity(mp);
    double kappa = BulkViscosity(mp);
    
    tens4ds c = dyad1s(I)*(kappa - 2.*mu/3.) + dyad4s(I)*(2*mu);
    
    return c;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to temperature T
mat3ds FENewtonianFluid::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    
    mat3ds D = vt.RateOfDeformation();
    
    double dmu = TangentShearViscosityTemperature(mp);
    double dkappa = TangentBulkViscosityTemperature(mp);
    
    mat3ds ds = mat3dd(1.0)*(D.tr()*(dkappa - 2.*dmu/3.)) + D*(2*dmu);
        
    return ds;
}

//-----------------------------------------------------------------------------
//! dynamic shear viscosity
double FENewtonianFluid::ShearViscosity(FEMaterialPoint& mp)
{
    double mu = m_mu;
    if (m_muhat) {
        FEThermoFluidMaterialPoint* tf = mp.ExtractData<FEThermoFluidMaterialPoint>();
        if (tf) {
            double That = (tf->m_T+m_Tr)/m_Tr;
            mu *= m_muhat->value(That);
        }
    }
    return mu;
}

//-----------------------------------------------------------------------------
//! dynamic shear viscosity tangent w.r.t. temperature
double FENewtonianFluid::TangentShearViscosityTemperature(FEMaterialPoint& mp)
{
    double dmu = 0;
    if (m_muhat) {
        FEThermoFluidMaterialPoint* tf = mp.ExtractData<FEThermoFluidMaterialPoint>();
        if (tf) {
            double That = (tf->m_T+m_Tr)/m_Tr;
            dmu = m_muhat->derive(That)*m_mu/m_Tr;
        }
    }
    return dmu;
}

//-----------------------------------------------------------------------------
//! bulk viscosity
double FENewtonianFluid::BulkViscosity(FEMaterialPoint& mp)
{
    double kappa = m_kappa;
    if (m_kappa) {
        FEThermoFluidMaterialPoint* tf = mp.ExtractData<FEThermoFluidMaterialPoint>();
        if (tf) {
            double That = (tf->m_T+m_Tr)/m_Tr;
            kappa *= m_kappahat->value(That);
        }
    }
    return kappa;
}

//-----------------------------------------------------------------------------
//! bulk viscosity tangent w.r.t. temperature
double FENewtonianFluid::TangentBulkViscosityTemperature(FEMaterialPoint& mp)
{
    double dkappa = 0;
    if (m_kappa) {
        FEThermoFluidMaterialPoint* tf = mp.ExtractData<FEThermoFluidMaterialPoint>();
        if (tf) {
            double That = (tf->m_T+m_Tr)/m_Tr;
            dkappa = m_kappahat->derive(That)*m_kappa/m_Tr;
        }
    }
    return dkappa;
}
