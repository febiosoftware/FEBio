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
#include "FENewtonianRealVapor.h"
#include "FEFluid.h"
#include "FEBiphasicFSI.h"
#include "FEThermoFluid.h"
#include "FERealVapor.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FENewtonianRealVapor, FEViscousFluid)
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"   )->setUnits(UNIT_VISCOSITY)->setLongName("referential shear viscosity");
// properties
    ADD_PROPERTY(m_esat , "esat",FEProperty::Optional)->SetLongName("saturation dilatation");
    ADD_PROPERTY(m_musat, "musat",FEProperty::Optional)->SetLongName("normalized saturation shear viscosity");
    ADD_PROPERTY(m_C[0] , "C0", FEProperty::Optional)->SetLongName("1st mu virial coeff");
    ADD_PROPERTY(m_C[1] , "C1", FEProperty::Optional)->SetLongName("2nd mu virial coeff");
    ADD_PROPERTY(m_C[2] , "C2", FEProperty::Optional)->SetLongName("3rd mu virial coeff");
    ADD_PROPERTY(m_C[3] , "C3", FEProperty::Optional)->SetLongName("4th mu virial coeff");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FENewtonianRealVapor::FENewtonianRealVapor(FEModel* pfem) : FEViscousFluid(pfem)
{
    m_kappa = 0;
    m_mu = 0;
    m_Tr = 0;
    m_nvc = 0;
    m_esat = nullptr;
    m_musat = nullptr;
    for (int k=0; k<MAX_NVC; ++k) m_C[k] = nullptr;
}

//-----------------------------------------------------------------------------
//! initialization
bool FENewtonianRealVapor::Init()
{
    if (m_musat) {
        m_Tr = GetGlobalConstant("T");
        
        if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    }
    
    if (m_esat) m_esat->Init();
    else {
        FECoreBase* pMat = GetAncestor();
        FEThermoFluid* pFluid = dynamic_cast<FEThermoFluid*>(pMat);
        if (pFluid) {
            FERealVapor* pRV = dynamic_cast<FERealVapor*>(pFluid->GetElastic());
            if (pRV) {
                m_esat = pRV->m_esat;
                m_esat->Init();
                m_Tc = pRV->m_Tc;
                m_alpha = pRV->m_alpha;
            }
            else return false;
        }
        else return false;
    }
    if (m_musat) m_musat->Init();
    m_nvc = 0;
    for (int k=0; k<MAX_NVC; ++k) {
        if (m_C[k]) {
            m_C[k]->Init();
            ++m_nvc;
        }
    }

    return FEViscousFluid::Init();
}

//-----------------------------------------------------------------------------
void FENewtonianRealVapor::Serialize(DumpStream& ar)
{
    FEViscousFluid::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    ar & m_Tr;
}

//-----------------------------------------------------------------------------
//! viscous stress
mat3ds FENewtonianRealVapor::Stress(FEMaterialPoint& pt)
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
mat3ds FENewtonianRealVapor::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = pf.RateOfDeformation();
    
    // evaluate ∂µ/∂J
    double dmuJ = TangentShearViscosityStrain(mp);
    double dkpJ = TangentBulkViscosityStrain(mp);

    mat3ds dsJ = mat3dd(1.0)*(D.tr()*(dkpJ - 2.*dmuJ/3.)) + D*(2*dmuJ);
    
    return dsJ;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to rate of deformation tensor D
tens4ds FENewtonianRealVapor::Tangent_RateOfDeformation(FEMaterialPoint& mp)
{
    mat3dd I(1.0);
    double mu = ShearViscosity(mp);
    double kappa = BulkViscosity(mp);
    
    tens4ds c = dyad1s(I)*(kappa - 2.*mu/3.) + dyad4s(I)*(2*mu);
    
    return c;
}

//-----------------------------------------------------------------------------
//! tangent of stress with respect to temperature T
mat3ds FENewtonianRealVapor::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();
    mat3ds D = pf.RateOfDeformation();

    // evaluate ∂µ/∂T
    double dmuT = TangentShearViscosityTemperature(mp);
    double dkpT = TangentBulkViscosityTemperature(mp);
    
    mat3ds dsT = mat3dd(1.0)*(D.tr()*(dkpT - 2.*dmuT/3.)) + D*(2*dmuT);
    
    return dsT;
}

//-----------------------------------------------------------------------------
//! dynamic shear viscosity
double FENewtonianRealVapor::ShearViscosity(FEMaterialPoint& mp)
{
    double mu = 1.0;
    if (m_musat) {
        FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
        FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();
        double T = tf.m_T + m_Tr;
        double That = T/m_Tr;
        double J = 1 + pf.m_ef;
        double y = (That < m_Tc) ? (m_Tc-That)/(m_Tc-1) : 0;
        double q = log(1+pow(y,m_alpha));
        mu = m_musat->value(q);
        double Jsat = exp(m_esat->value(q));
        double x = 1 - Jsat/J;
        for (int k=0; k<m_nvc; ++k) mu += m_C[k]->value(q)*pow(x,k+1);
    }
    return mu*m_mu;
}

//-----------------------------------------------------------------------------
//! dynamic shear viscosity tangent w.r.t. temperature
double FENewtonianRealVapor::TangentShearViscosityTemperature(FEMaterialPoint& mp)
{
    double d = 1e-6*m_Tr;
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();
    
    // evaluate ∂µ/∂T
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = pf.m_ef;
    ft->m_T = tf.m_T+d;
    FEMaterialPoint tmp(ft);
    double mup = ShearViscosity(tmp);
    ft->m_T = tf.m_T-d;
    double mum = ShearViscosity(tmp);
    delete ft;
    double dmuT = (mup - mum)/(2*d);
    
    return dmuT;
}

//-----------------------------------------------------------------------------
//! dynamic shear viscosity tangent w.r.t. J
double FENewtonianRealVapor::TangentShearViscosityStrain(FEMaterialPoint& mp)
{
    double d = 1e-6;
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();
    
    // evaluate ∂µ/∂J
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = pf.m_ef+d;
    ft->m_T = tf.m_T;
    FEMaterialPoint tmp(ft);
    double mup = ShearViscosity(tmp);
    fp->m_ef = pf.m_ef-d;
    double mum = ShearViscosity(tmp);
    delete ft;
    double dmuJ = (mup - mum)/(2*d);

    return dmuJ;
}

//-----------------------------------------------------------------------------
//! bulk viscosity
double FENewtonianRealVapor::BulkViscosity(FEMaterialPoint& mp)
{
    double kappa = 0;
    return kappa;
}

//-----------------------------------------------------------------------------
//! bulk viscosity tangent w.r.t. temperature
double FENewtonianRealVapor::TangentBulkViscosityTemperature(FEMaterialPoint& mp)
{
    double dkpT = 0;
    return dkpT;
}

//-----------------------------------------------------------------------------
//! bulk viscosity tangent w.r.t. J
double FENewtonianRealVapor::TangentBulkViscosityStrain(FEMaterialPoint& mp)
{
    double dkpJ = 0;
    return dkpJ;
}
