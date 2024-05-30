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
ADD_PROPERTY(m_D[0] , "D0"   )->SetLongName("1st normalized pressure coefficient");
ADD_PROPERTY(m_D[1] , "D1" , FEProperty::Optional)->SetLongName("2nd normalized pressure coefficient");
ADD_PROPERTY(m_D[2] , "D2" , FEProperty::Optional)->SetLongName("3rd normalized pressure coefficient");
ADD_PROPERTY(m_D[3] , "D3" , FEProperty::Optional)->SetLongName("4th normalized pressure coefficient");
ADD_PROPERTY(m_D[4] , "D4" , FEProperty::Optional)->SetLongName("5th normalized pressure coefficient");
ADD_PROPERTY(m_D[5] , "D5" , FEProperty::Optional)->SetLongName("6th normalized pressure coefficient");
ADD_PROPERTY(m_D[6] , "D6" , FEProperty::Optional)->SetLongName("7th normalized pressure coefficient");

ADD_PROPERTY(m_C[0], "C0")->SetLongName("1st cv virial coeff");
ADD_PROPERTY(m_C[1], "C1", FEProperty::Optional)->SetLongName("2nd cv virial coeff");
ADD_PROPERTY(m_C[2], "C2", FEProperty::Optional)->SetLongName("3rd cv virial coeff");
ADD_PROPERTY(m_C[3], "C3", FEProperty::Optional)->SetLongName("4th cv virial coeff");

END_FECORE_CLASS();

FERealVapor::FERealVapor(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_nvp = 0;
    m_nvc = 0;
    m_Pr = m_Tr = 0;
    for (int k=0; k<MAX_NVP; ++k) m_D[k] = nullptr;
    for (int k=0; k<MAX_NVC; ++k) m_C[k] = nullptr;
    m_psat = m_asat = m_ssat = m_esat = m_cvsat = nullptr;
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
        if (m_D[k]) {
            m_D[k]->Init();
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
//! gauge pressure
double FERealVapor::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double J = 1 + fp.m_ef;
    double Jsat = 1 + m_esat->value(That);
    double Psat = 1 + m_psat->value(That);
    double D[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) D[k] = m_D[k]->value(That);
    double x = Jsat/J;
    double sum = 0;
    for (int k=0; k<m_nvp; ++k) sum += D[k]*pow(x,k+1);
    double p = Psat*(x + (1-x)*sum) - 1;
    
    return p*m_Pr;
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
    double D[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) D[k] = m_D[k]->value(That);
    double x = Jsat/J;
    double sum = 0;
    for (int k=1; k<m_nvp; ++k) {
        double xk = pow(x,k);
        double xkp = xk*x;
        sum += D[k]*((1-xkp)/(k+1) - (1-xk)/k);
    }
    double a = asat - Jsat*(1-1/x) + Psat*Jsat*(log(x) + D[0]*(log(x)+1-x) + sum);

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
    
    double J = 1 + fp.m_ef;
    double Jsat = 1 + m_esat->value(That);
    double Psat = 1 + m_psat->value(That);
    double dJsat = m_esat->derive(That);
    double dPsat = m_psat->derive(That);
    double D[MAX_NVP], dD[MAX_NVP];
    for (int k=0; k<m_nvp; ++k) {
        D[k] = m_D[k]->value(That);
        dD[k] = m_D[k]->derive(That);
    }
    double x = Jsat/J;
    double dx = dJsat/J;
    double sum1 = 0, sum2 = 0, sum3 = 0;
    for (int k=1; k<m_nvp; ++k) {
        double xk = pow(x,k);
        double xkp = xk*x;
        double y = (1-xkp)/(k+1) - (1-xk)/k;
        sum1 += D[k]*y;
        sum2 += D[k]*xk;
        sum3 += dD[k]*y;
    }
    
    double z = log(x)+1-x;
    double s = ssat - (dPsat*Jsat + Psat*dJsat)*(log(x) + D[0]*z + sum1)
    - Psat*(1-x)*dJsat*(D[0] + sum2) - Psat*Jsat*(dD[0]*z + sum3);

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
    for (int k=0; k<m_nvc; ++k) cv += m_C[k]->value(That)*pow(y,k+1);

    return cv*m_Pr/(m_Tr*m_rhor);
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
    double T = m_Tr + tf.m_T;
    double cp = cv - T/m_rhor*pow(dpT,2)/dpJ;
    return cp;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FERealVapor::Dilatation(const double T, const double p, double& e)
{
    // check that the prescribed gauge pressure is valid (positive)
    if (p < 0) return false;
    // if valid, continue
    double Phat = 1 + p/m_Pr;
    double That = (T+m_Tr)/m_Tr;
    double Psat = 1 + m_psat->value(That);
    // check to make sure that we are in the vapor phase
    if (Phat > Psat) return false;
    // then continue
    double Jsat = 1 + m_esat->value(That);
    double D[MAX_NVP];
    vector <double> B(m_nvp+2,0);
    for (int k=0; k<m_nvp; ++k) D[k] = m_D[k]->value(That);
    B[0] = -Phat/Psat;
    B[1] = 1 + D[0];
    B[m_nvp+1] = -D[m_nvp-1];
    for (int k=0; k<m_nvp-1; ++k) B[k+2] = D[k+1] - D[k];
    double x = (-B[1]+sqrt(B[1]*B[1]-4*B[0]*B[2]))/(2*B[2]);
    // check to see if we are infinitesimally close to the saturation curve
    if (Psat - Phat < 1e-3) {
        // if nearly on saturation curve, use saturation curve dilatation
        e = Jsat/x - 1;
        return true;
    }
    // initial guess for J depends if e = 0 or not
    double J = (e == 0) ? Jsat/x : 1 + e;
    // if one virial coefficient only, we're done
    if (m_nvp == 1) { e = J - 1; return true; }
    // solve iteratively for J using Newton's method
    bool convgd = solvepoly(m_nvp, B, x, true);
    J = Jsat/x;
    e = J - 1;
    return convgd;
}
