/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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


#include "FERealLiquid.h"
#include "FEFluidMaterialPoint.h"
#include "FEThermoFluidMaterialPoint.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERealLiquid, FEElasticFluid)

// material parameters
ADD_PARAMETER(m_nvc  , FE_RANGE_GREATER(0), "nvc");

ADD_PROPERTY(m_psat, "psat")->SetLongName("saturation gauge pressure normalized");
ADD_PROPERTY(m_asat, "asat")->SetLongName("saturation free energy normalized");
ADD_PROPERTY(m_ssat, "ssat")->SetLongName("saturation entropy normalized");
ADD_PROPERTY(m_esat, "esat")->SetLongName("saturation dilatation");
ADD_PROPERTY(m_B[0], "B1"  )->SetLongName("1st pressure virial coefficient");
ADD_PROPERTY(m_B[1], "B2", FEProperty::Optional)->SetLongName("2nd pressure virial coefficient");
ADD_PROPERTY(m_B[2], "B3", FEProperty::Optional)->SetLongName("3rd pressure virial coefficient");
ADD_PROPERTY(m_cvsat, "cvsat")->SetLongName("saturation isochoric heat capacity normalized");
ADD_PROPERTY(m_C[0], "C1"  )->SetLongName("1st cv virial coefficient");
ADD_PROPERTY(m_C[1], "C2", FEProperty::Optional)->SetLongName("2nd cv virial coefficient");
ADD_PROPERTY(m_C[2], "C3", FEProperty::Optional)->SetLongName("3rd cv virial coefficient");

END_FECORE_CLASS();

FERealLiquid::FERealLiquid(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_nvc = 0;
    m_R = m_Pr = m_Tr = 0;
    m_psat = m_asat = m_ssat = m_esat = m_cvsat = nullptr;
    m_B[0] = m_B[1] = m_B[2] = nullptr;
    m_C[0] = m_C[1] = m_C[2] = nullptr;
}

//-----------------------------------------------------------------------------
//! initialization
bool FERealLiquid::Init()
{
    m_R  = GetGlobalConstant("R");
    m_Tr = GetGlobalConstant("T");
    m_Pr = GetGlobalConstant("P");
    
    if (m_R  <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section");    return false; }
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    if (m_Pr <= 0) { feLogError("A positive referential absolute pressure P must be defined in Globals section"); return false; }
    if (m_nvc == 0){ feLogError("At least one virial coefficient must be specified in this real liquid"); return false; }

    m_pMat = dynamic_cast<FEThermoFluid*>(GetParent());
    m_rhor = m_pMat->ReferentialDensity();
    
    m_psat->Init();
    m_asat->Init();
    m_ssat->Init();
    m_esat->Init();
    m_cvsat->Init();
    for (int k=0; k<m_nvc; ++k) {
        if (m_B[k]) m_B[k]->Init();
        if (m_C[k]) m_C[k]->Init();
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FERealLiquid::Serialize(DumpStream& ar)
{
    FEElasticFluid::Serialize(ar);
    
    if (ar.IsShallow()) return;
    ar & m_pMat;
    ar & m_R & m_Pr & m_Tr & m_rhor;
}

//-----------------------------------------------------------------------------
//! gage pressure
double FERealLiquid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double p = m_psat->value(That);
    double x = fp.m_ef - m_esat->value(That);
    for (int k=1; k<=m_nvc; ++k)
        p += m_B[k-1]->value(That)*pow(x,k);
    
    return p*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FERealLiquid::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double dpJ = m_B[0]->value(That);
    double x = fp.m_ef - m_esat->value(That);
    for (int k=2; k<=m_nvc; ++k)
        dpJ += k*m_B[k-1]->value(That)*pow(x,k-1);
    
    return dpJ*m_Pr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FERealLiquid::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double dpJ2 = (m_nvc > 1) ? 2*m_B[1]->value(That) : 0;
    double x = fp.m_ef - m_esat->value(That);
    for (int k=3; k<=m_nvc; ++k)
        dpJ2 += k*(k-1)*m_B[k-1]->value(That)*pow(x,k-2);
    
    return dpJ2*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FERealLiquid::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double dpT = m_psat->derive(That);
    double dJsat = m_esat->derive(That);
    double x = fp.m_ef - m_esat->value(That);
    for (int k=1; k<=m_nvc; ++k)
        dpT += m_B[k-1]->derive(That)*pow(x,k)
        - k*dJsat*m_B[k-1]->value(That)*pow(x,k-1);
    
    return dpT*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FERealLiquid::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double dpT2 = m_psat->deriv2(That);
    double dJsat = m_esat->derive(That);
    double dJsa2 = m_esat->deriv2(That);
    double x = fp.m_ef - m_esat->value(That);
    vector<double> B(m_nvc,0);
    for (int k=0; k<m_nvc; ++k) B[k] = m_B[k]->value(That);
    for (int k=1; k<=m_nvc; ++k)
        dpT2 += m_B[k-1]->deriv2(That)*pow(x,k)
        - (2*dJsat*m_B[k-1]->derive(That) + dJsa2*B[k-1])*k*pow(x,k-1);
    for (int k=2; k<=m_nvc; ++k)
        dpT2 += k*(k-1)*pow(dJsat,2)*B[k-1]*pow(x,k-2);

    return dpT2*m_Pr/pow(m_Tr,2);
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FERealLiquid::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double dpJT = m_B[0]->derive(That);
    double dJsat = m_esat->derive(That);
    double x = fp.m_ef - m_esat->value(That);
    for (int k=2; k<=m_nvc; ++k)
        dpJT += k*m_B[k-1]->derive(That)*pow(x,k-1)
        - k*(k-1)*dJsat*m_B[k-1]->value(That)*pow(x,k-2);
    
    return dpJT*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FERealLiquid::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double x = fp.m_ef - m_esat->value(That);
    double a = m_asat->value(That) - x*m_psat->value(That);
    for (int k=1; k<=m_nvc; ++k)
        a -= m_B[k-1]->value(That)/(k+1)*pow(x,k+1);
    
    return a*m_Pr/m_rhor;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FERealLiquid::SpecificEntropy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double x = fp.m_ef - m_esat->value(That);
    double dJsat = m_esat->derive(That);
    double ssat = m_ssat->value(That);
    double s = ssat + x*m_psat->derive(That);
    for (int k=1; k<=m_nvc; ++k)
        s += (m_B[k-1]->derive(That)*x/(k+1) - dJsat*m_B[k-1]->value(That))*pow(x,k);

    return s*m_Pr/(m_rhor*m_Tr);
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FERealLiquid::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();

    // get the specific free energy
    double a = SpecificFreeEnergy(mp);
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double asat = m_asat->value(That)*m_Pr/m_rhor;
    
    // the specific strain energy is the difference between these two values
    return a - asat;
}

//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FERealLiquid::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double x = fp.m_ef - m_esat->value(That);
    double cv = 0;
    for (int k=1; k<=m_nvc; ++k) cv += m_C[k-1]->value(That)*pow(x,k);
    cv *= m_cvsat->value(1.0);
    cv += m_cvsat->value(That);
    
    return m_Pr/(m_rhor*m_Tr)*cv;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FERealLiquid::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double x = fp.m_ef - m_esat->value(That);
    double dcvJ = 0;
    for (int k=1; k<=m_nvc; ++k) dcvJ += k*m_C[k-1]->value(That)*pow(x,k-1);
    dcvJ *= m_cvsat->value(1.0);
    
    return m_Pr/(m_rhor*m_Tr)*dcvJ;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FERealLiquid::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double That = (m_Tr+tf.m_T)/m_Tr;
    double x = fp.m_ef - m_esat->value(That);
    double dcvT = 0;
    for (int k=1; k<=m_nvc; ++k) dcvT += (m_C[k-1]->derive(That)*x - k*m_esat->derive(That)*m_C[k-1]->value(That))*pow(x,k-1);
    dcvT *= m_cvsat->value(1.0);
    dcvT += m_cvsat->derive(That);

    return m_Pr/(m_rhor*pow(m_Tr,2))*dcvT;
}

//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FERealLiquid::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double cv = IsochoricSpecificHeatCapacity(mp);
    double p = Pressure(mp);
    // temporarily reduce the number of virial coefficients to 1 for the next calculation
    int nvc = m_nvc;
    m_nvc = 1;
    double dpT = Tangent_Temperature(mp);
    double dpJ = Tangent_Strain(mp);
    m_nvc = nvc;
    double cp = cv + (dpT/dpJ)/m_rhor*(p - (m_Tr + tf.m_T)*dpT);
    return cp;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FERealLiquid::Dilatation(const double T, const double p, double& e)
{
    double errrel = 1e-6;
    double errabs = 1e-6;
    int maxiter = 100;
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    double That = (T+m_Tr)/m_Tr;
    double psat = m_psat->value(That);
    double esat = m_esat->value(That);
    double phat = p/m_Pr;
    if (phat == psat) { e = esat; return true; }
    switch (m_nvc) {
        case 0:
            delete ft;
            return false;
            break;
        case 1:
        {
            double B0 = m_B[0]->value(That);
            if (B0 == 0) { e = esat; return true; }
            e = (phat - psat)/B0 + esat;
            delete ft;
            return true;
        }
            break;
        case 2:
        {
            double B0 = m_B[0]->value(That);
            double B1 = m_B[1]->value(That);
            if ((B0 == 0) || (B1 == 0)) { e = esat; return true; }
            double det = pow(B0,2) + 4*B1*(phat - psat);
            if (det < 0) return false;
            det = sqrt(det);
            // only one root of this quadratic equation is valid.
            e = esat - (B0+det)/(2*B1);
            delete ft;
            return true;
        }
            break;

        default:
        {
            bool convgd = false;
            bool done = false;
            int iter = 0;
			FEMaterialPoint mp(ft);
			do {
                ++iter;
                fp->m_ef = e;
                double f = Pressure(mp) - p;
                double df = Tangent_Strain(mp);
                double de = (df != 0) ? -f/df : 0;
                e += de;
                if ((fabs(de) < errrel*fabs(e)) ||
                    (fabs(f/m_Pr) < errabs)) { convgd = true; done = true; }
                if (iter > maxiter) done = true;
            } while (!done);
            
            delete ft;
            return convgd;
        }
            break;
    }
    delete ft;
    return false;
}

//-----------------------------------------------------------------------------
//! pressure from state variables
double FERealLiquid::Pressure(const double ef, const double T)
{
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = ef;
    ft->m_T = T;
    FEMaterialPoint tmp(ft);
    double p = Pressure(tmp);
    delete ft;
    return p;
}
