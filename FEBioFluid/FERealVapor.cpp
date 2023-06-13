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


#include "FERealVapor.h"
#include <FECore/log.h>
#include <FECore/FEFunction1D.h>
#include <FECore/sys.h>
#include "FEFluidMaterialPoint.h"
#include "FEThermoFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERealVapor, FEElasticFluid)

// material parameters
ADD_PROPERTY(m_esat , "esat" )->SetLongName("saturation dilatation");
ADD_PROPERTY(m_psat , "psat" )->SetLongName("saturation gauge pressure normalized");
ADD_PROPERTY(m_asat , "asat" )->SetLongName("saturation free energy normalized");
ADD_PROPERTY(m_ssat , "ssat" )->SetLongName("saturation entropy normalized");
ADD_PROPERTY(m_cvsat, "cvsat")->SetLongName("saturation cv normalized");
ADD_PROPERTY(m_B[0] , "B1"   )->SetLongName("1st normalized pressure coefficient");
ADD_PROPERTY(m_B[1] , "B2" , FEProperty::Optional)->SetLongName("2nd normalized pressure coefficient");
ADD_PROPERTY(m_B[2] , "B3" , FEProperty::Optional)->SetLongName("3rd normalized pressure coefficient");
ADD_PROPERTY(m_B[3] , "B4" , FEProperty::Optional)->SetLongName("4th normalized pressure coefficient");
ADD_PROPERTY(m_B[4] , "B5" , FEProperty::Optional)->SetLongName("5th normalized pressure coefficient");
ADD_PROPERTY(m_B[5] , "B6" , FEProperty::Optional)->SetLongName("6th normalized pressure coefficient");
ADD_PROPERTY(m_B[6] , "B7" , FEProperty::Optional)->SetLongName("7th normalized pressure coefficient");

ADD_PROPERTY(m_C[0], "C0")->SetLongName("1st cv virial coeff");
ADD_PROPERTY(m_C[1], "C1", FEProperty::Optional)->SetLongName("2nd cv virial coeff");
ADD_PROPERTY(m_C[2], "C2", FEProperty::Optional)->SetLongName("3rd cv virial coeff");

END_FECORE_CLASS();

FERealVapor::FERealVapor(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_nvp = 0;
    m_nvc = 0;
    m_Pr = m_Tr = 0;
    for (int k=0; k<MAX_NVP; ++k) m_B[k] = nullptr;
    for (int k=0; k<MAX_NVC; ++k) m_C[k] = nullptr;
    m_psat = m_asat = m_ssat = nullptr;
}

//-----------------------------------------------------------------------------
//! initialization
bool FERealVapor::Init()
{
    m_Tr = GetGlobalConstant("T");
    m_Pr = GetGlobalConstant("P");
    
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    if (m_Pr <= 0) { feLogError("A positive referential absolute pressure P must be defined in Globals section"); return false; }

    m_pMat = dynamic_cast<FEThermoFluid*>(GetParent());
    m_rhor = m_pMat->ReferentialDensity();

    m_esat->Init();
    m_psat->Init();
    m_asat->Init();
    m_ssat->Init();
    m_cvsat->Init();
    m_nvp = 0;
    for (int k=0; k<MAX_NVP; ++k) {
        if (m_B[k]) {
            m_B[k]->Init();
            ++m_nvp;
        }
    }
    m_nvc = 0;
    for (int k=0; k<MAX_NVC; ++k) {
        if (m_C[k]) {
            m_C[k]->Init();
            ++m_nvc;
        }
    }

    if (m_nvp < 1) { feLogError("At least one virial coefficient should be provided for the pressure"); return false; }
    if (m_nvc < 1) { feLogError("At least one virial coefficient should be provided for cv"); return false; }
    
    return true;
}

//-----------------------------------------------------------------------------
void FERealVapor::Serialize(DumpStream& ar)
{
    FEElasticFluid::Serialize(ar);

    if (ar.IsShallow()) return;
    ar & m_pMat;
    ar & m_Pr & m_Tr & m_rhor;
    ar & m_nvp & m_nvc;
}

//-----------------------------------------------------------------------------
//! gage pressure
double FERealVapor::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double J = 1 + fp.m_ef;
    double Jsat = 1 + m_esat->value(That);
    double Psat = 1 + m_psat->value(That);
    double B[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) B[k] = m_B[k]->value(That);
    double x = log(J/Jsat);
    double y = 0;
    for (int k=0; k<m_nvp; ++k) y += B[k]*pow(x,k+1);
    double P = Psat*exp(y);
    
    return (P-1)*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FERealVapor::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double J = 1 + fp.m_ef;
    double Jsat = 1 + m_esat->value(That);
    double Psat = 1 + m_psat->value(That);
    double B[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) B[k] = m_B[k]->value(That);
    double x = log(J/Jsat);
    double y = 0;
    for (int k=0; k<m_nvp; ++k) y += B[k]*pow(x,k+1);
    double P = Psat*exp(y);
    double df = 0;
    for (int k=0; k<m_nvp; ++k) df += (k+1)*B[k]*pow(x,k);
    df /= J;
    double dpJ = P*df;

    return dpJ*m_Pr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FERealVapor::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double J = 1 + fp.m_ef;
    double Jsat = 1 + m_esat->value(That);
    double Psat = 1 + m_psat->value(That);
    double B[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) B[k] = m_B[k]->value(That);
    double x = log(J/Jsat);
    double y = 0;
    for (int k=0; k<m_nvp; ++k) y += B[k]*pow(x,k+1);
    double P = Psat*exp(y);
    double df = 0, d2f = -B[0];
    for (int k=0; k<m_nvp; ++k) df += (k+1)*B[k]*pow(x,k);
    for (int k=1; k<m_nvp; ++k) d2f += (k+1)*B[k]*pow(x,k-1)*(k - x);
    df /= J;
    d2f /= J*J;
    double dpJ2 = P*(d2f + df*df);

    return dpJ2*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FERealVapor::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double J = 1 + fp.m_ef;
    double Jsat = 1 + m_esat->value(That);
    double dJsat = m_esat->derive(That);
    double Psat = 1 + m_psat->value(That);
    double dPsat = m_psat->derive(That);
    double B[MAX_NVP], dB[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) {
        B[k] = m_B[k]->value(That);
        dB[k] = m_B[k]->derive(That);
    }
    double x = log(J/Jsat);
    double y = 0;
    for (int k=0; k<m_nvp; ++k) y += B[k]*pow(x,k+1);
    double P = Psat*exp(y);
    double df = dPsat/P;
    for (int k=0; k<m_nvp; ++k) df += pow(x,k+1)*(dB[k] - (k+1)*B[k]*dJsat/(Jsat*x));
    double dpT = P*df;
    
    return dpT*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FERealVapor::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double J = 1 + fp.m_ef;
    double Jsat = 1 + m_esat->value(That);
    double dJsat = m_esat->derive(That);
    double d2Jsat = m_esat->deriv2(That);
    double Psat = 1 + m_psat->value(That);
    double dPsat = m_psat->derive(That);
    double d2Psat = m_psat->deriv2(That);
    double B[MAX_NVP], dB[MAX_NVP], d2B[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) {
        B[k] = m_B[k]->value(That);
        dB[k] = m_B[k]->derive(That);
        d2B[k] = m_B[k]->deriv2(That);
    }
    double x = log(J/Jsat);
    double y = 0;
    for (int k=0; k<m_nvp; ++k) y += B[k]*pow(x,k+1);
    double P = Psat*exp(y);
    double df = dPsat/P;
    double d2f = (d2Psat - dPsat*df)/P;
    for (int k=0; k<m_nvp; ++k) {
        double px = pow(x,k+1);
        df += px*(dB[k] - (k+1)*B[k]*dJsat/(Jsat*x));
        d2f += px*d2B[k] - (k+1)*pow(x,k)*(2*dB[k]*dJsat + B[k]*d2Jsat)/Jsat;
    }
    for (int k=1; k<m_nvp; ++k) d2f += (k+1)*pow(x,k-1)*(k+x)*B[k]*pow(dJsat/Jsat,2);
    double dpT2 = P*(d2f + df*df);

    return dpT2*m_Pr/pow(m_Tr,2);
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FERealVapor::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double J = 1 + fp.m_ef;
    double Jsat = 1 + m_esat->value(That);
    double Psat = 1 + m_psat->value(That);
    double dJsat = m_esat->derive(That);
    double dPsat = m_psat->derive(That);
    double B[MAX_NVP], dB[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) {
        B[k] = m_B[k]->value(That);
        dB[k] = m_B[k]->derive(That);
    }
    double x = log(J/Jsat);
    double y = 0;
    for (int k=0; k<m_nvp; ++k) y += B[k]*pow(x,k+1);
    double P = Psat*exp(y);
    double dfJ = 0;
    double dfT = dPsat/P;
    double dfTJ = dB[0];
    for (int k=0; k<m_nvp; ++k) {
        dfJ += (k+1)*B[k]*pow(x,k);
        dfT += pow(x,k+1)*(dB[k] - (k+1)*B[k]*dJsat/(Jsat*x));
    }
    dfJ /= J;
    for (int k=1; k<m_nvp; ++k) dfTJ += (k+1)*pow(x,k-1)*(dB[k]*x - k*B[k]*dJsat/Jsat);
    dfTJ /= J;
    double dpTJ = P*(dfTJ + dfT*dfJ);

    return dpTJ*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FERealVapor::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double asat = m_asat->value(That);
    // incomplete
    double a = asat;

    double J = 1 + fp.m_ef;
    double Jsat = 1 + m_esat->value(That);
    double Psat = 1 + m_psat->value(That);
    double B[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) B[k] = m_B[k]->value(That);
    double dJ = 1e-3;
    double j = Jsat;
    double p_p = Psat - 1, p_t = 0;
    bool done = false, penult = false;
    do {
        if (penult) done = true;
        double x = log(j/Jsat);
        double y = 0;
        for (int k=0; k<m_nvp; ++k) y += B[k]*pow(x,k+1);
        p_t = Psat*exp(y) - 1;
        a -= 0.5*dJ*(p_p + p_t);
        if (j + dJ > J) { dJ = J - j; penult = true; }
        j += dJ;
        p_p = p_t;
    } while (!done);

    return a*m_Pr/m_rhor;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FERealVapor::SpecificEntropy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double ssat = m_ssat->value(That);
    // incomplete
    double s = ssat;

    return s*m_Pr/(m_Tr*m_rhor);
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FERealVapor::SpecificStrainEnergy(FEMaterialPoint& mp)
{
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
double FERealVapor::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();

    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double Jsat = 1 + m_esat->value(That);
    double y = 1 - Jsat/J;
    double cv = m_cvsat->value(That);
    for (int k=0; k<m_nvc; ++k) cv -= m_C[k]->value(That)*pow(y,k);

    return cv*m_Pr/(m_Tr*m_rhor);
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FERealVapor::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double Jsat = 1 + m_esat->value(That);
    double y = 1 - Jsat/J;
    double dcvJ = 0;
    for (int k=0; k<m_nvc; ++k) dcvJ -= (k+1)*m_C[k]->value(That)*pow(y,k);
    dcvJ *= Jsat/(J*J);
    
    return dcvJ*m_Pr/(m_Tr*m_rhor);
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FERealVapor::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double Jsat = 1 + m_esat->value(That);
    double dJsat = m_esat->derive(That);
    double y = 1 - Jsat/J;
    double dcvT = m_cvsat->derive(That);
    for (int k=0; k<m_nvc; ++k) dcvT -= (m_C[k]->derive(That)*y - (k+1)*dJsat/J*m_C[k]->value(That))*pow(y,k);
    return dcvT*m_Pr/(pow(m_Tr,2)*m_rhor);
}

//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FERealVapor::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double cv = IsochoricSpecificHeatCapacity(mp);
    double p = Pressure(mp);
    double dpT = Tangent_Temperature(mp);
    double dpJ = Tangent_Strain(mp);
    double cp = cv + (dpT/dpJ)/m_rhor*(p - (m_Tr + tf.m_T)*dpT);
    return cp;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FERealVapor::Dilatation(const double T, const double p, double& e)
{
    double Phat = 1 + p/m_Pr;
    double That = (T+m_Tr)/m_Tr;
    double Psat = 1 + m_psat->value(That);
    double Jsat = 1 + m_esat->value(That);
    double B[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) B[k] = m_B[k]->value(That);
    // initial guess for J depends if e = 0 or not
    double J = (e == 0) ? Jsat*pow(Phat/Psat,1/B[0]) : 1 + e;
    // if one virial coefficient only, we're done
    if (m_nvc == 1) { e = J - 1; return true; }
    // solve iteratively for J using Newton's method
    double errabs = 1e-12, errrel = 1e-6;
    int iter = 0, maxiter = 1000;
    double x = log(J/Jsat);
    bool convgd = false;
    do {
        double f = log(Phat/Psat);
        double df = 0;
        for (int k=0; k<m_nvc; ++k) {
            double pxk = pow(x,k);
            f -= B[k]*pxk*x;
            df -= (k+1)*B[k]*pxk;
        }
        double dx = (fabs(df) > 0) ? -f/df : 0;
        x += dx;
        if (fabs(dx) <= errrel*fabs(x)) convgd = true;
        if (fabs(f) <= errabs) convgd = true;
        ++iter;
    } while ((!convgd) && (iter < maxiter));
    J = Jsat*exp(x);
    e = J - 1;
    return convgd;
}
