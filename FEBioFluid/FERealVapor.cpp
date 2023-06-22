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
#include <FECore/tools.h>
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
ADD_PROPERTY(m_B[0] , "B0"   )->SetLongName("0th normalized pressure coefficient");
ADD_PROPERTY(m_B[1] , "B1" , FEProperty::Optional)->SetLongName("1st normalized pressure coefficient");
ADD_PROPERTY(m_B[2] , "B2" , FEProperty::Optional)->SetLongName("2nd normalized pressure coefficient");
ADD_PROPERTY(m_B[3] , "B3" , FEProperty::Optional)->SetLongName("3rd normalized pressure coefficient");
ADD_PROPERTY(m_B[4] , "B4" , FEProperty::Optional)->SetLongName("4th normalized pressure coefficient");
ADD_PROPERTY(m_B[5] , "B5" , FEProperty::Optional)->SetLongName("5th normalized pressure coefficient");
ADD_PROPERTY(m_B[6] , "B6" , FEProperty::Optional)->SetLongName("6th normalized pressure coefficient");

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
    double x = Jsat/J;
    double y = 0;
    for (int k=0; k<m_nvp; ++k) y += B[k]*pow(x,k+1);
    double P = Psat*y;
    
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
    double x = Jsat/J;
    double y = 0;
    for (int k=0; k<m_nvp; ++k) y += B[k]*(k+1)*pow(x,k+1);
    double dpJ = -Psat*y/J;

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
    double x = Jsat/J;
    double y = 0;
    for (int k=0; k<m_nvp; ++k) y += B[k]*(k+1)*(k+2)*pow(x,k+1);
    double dpJ2 = Psat*y/pow(J,2);

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
    double x = Jsat/J;
    double dx = dJsat/J;
    double dpT = 0;
    for (int k=0; k<m_nvp; ++k)
        dpT += (dPsat*B[k]+Psat*dB[k])*pow(x,k+1)
        + Psat*(k+1)*B[k]*dx*pow(dx,k);
    
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
    double x = Jsat/J;
    double dx = dJsat/J;
    double d2x = d2Jsat/J;
    double dpT2 = 0;
    for (int k=0; k<m_nvp; ++k) {
        dpT2 += (d2Psat*B[k]+2*dPsat*dB[k]+Psat*d2B[k])*pow(x,k+1)
        + (k+1)*(2*(dPsat*B[k]+Psat*dB[k])*dx
                 +Psat*B[k]*d2x)*pow(x,k);
    }
    for (int k=1; k<m_nvp; ++k)
        dpT2 += k*(k+1)*Psat*B[k]*pow(dx,2)*pow(x,k-1);

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
    double x = Jsat/J;
    double dx = dJsat/J;
    double dpTJ = 0;
    for (int k=0; k<m_nvp; ++k) {
        dpTJ -= (k+1)*pow(x,k+1)*(dPsat*B[k]+Psat*dB[k]
                                  +(k+1)*Psat*B[k]*dx);
    }
    dpTJ /= J;

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

    double J = 1 + fp.m_ef;
    double Jsat = 1 + m_esat->value(That);
    double Psat = 1 + m_psat->value(That);
    double B[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) B[k] = m_B[k]->value(That);
    double x = Jsat/J;
    double a = B[0]*log(x);
    for (int k=1; k<m_nvp; ++k)
        a -= B[k]/k*(1-pow(x,k));
    a = Psat*Jsat*a + asat + J - Jsat;

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
    double x = Jsat/J;
    double dx = dJsat/J;
    double y = dPsat/Psat + dJsat/Jsat;
    double s = m_ssat->value(That) + dJsat*Psat*(1-B[0]) - Jsat*Psat*(y*B[0]+dB[0])*log(x);
    for (int k=1; k<m_nvp; ++k) {
        double px = pow(x,k);
        s += Psat*(-dJsat*B[k]*px + Jsat*(y*B[k]+dB[k])/k*(1-px));
    }

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
    vector <double> B(m_nvp+1,0);
    for (int k=0; k<m_nvp; ++k) B[k+1] = m_B[k]->value(That);
    B[0] = -Phat/Psat;
    // initial guess for J depends if e = 0 or not
    double J = (e == 0) ? B[1]*Jsat*Psat/Phat : 1 + e;
    // if one virial coefficient only, we're done
    if (m_nvp == 1) { e = J - 1; return true; }
    // solve iteratively for J using Newton's method
    double x = Jsat/J;
    bool convgd = solvepoly(m_nvp, B, x, false);
    J = Jsat/x;
    e = J - 1;
    return convgd;
}
