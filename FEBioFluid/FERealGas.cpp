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


#include "FERealGas.h"
#include <FECore/log.h>
#include <FECore/FEFunction1D.h>
#include "FEFluidMaterialPoint.h"
#include "FEThermoFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERealGas, FEElasticFluid)

// material parameters
ADD_PARAMETER(m_nva  , FE_RANGE_GREATER_OR_EQUAL(0), "nva")->setLongName("no. of p virial coeffs");
ADD_PARAMETER(m_nvc  , FE_RANGE_GREATER_OR_EQUAL(0), "nvc")->setLongName("no. of cv virial coeffs");

ADD_PROPERTY(m_a0  , "a0")->SetLongName("normalized specific free energy 0");
ADD_PROPERTY(m_A[0], "A1", FEProperty::Optional)->SetLongName("1st p virial coeff");
ADD_PROPERTY(m_A[1], "A2", FEProperty::Optional)->SetLongName("2nd p virial coeff");
ADD_PROPERTY(m_A[2], "A3", FEProperty::Optional)->SetLongName("3rd p virial coeff");
ADD_PROPERTY(m_A[3], "A4", FEProperty::Optional)->SetLongName("4th p virial coeff");
ADD_PROPERTY(m_A[4], "A5", FEProperty::Optional)->SetLongName("5th p virial coeff");
ADD_PROPERTY(m_A[5], "A6", FEProperty::Optional)->SetLongName("6th p virial coeff");
ADD_PROPERTY(m_A[6], "A7", FEProperty::Optional)->SetLongName("7th p virial coeff");

ADD_PROPERTY(m_C[0], "C0")->SetLongName("1st cv virial coeff");
ADD_PROPERTY(m_C[1], "C1", FEProperty::Optional)->SetLongName("2nd cv virial coeff");
ADD_PROPERTY(m_C[2], "C2", FEProperty::Optional)->SetLongName("3rd cv virial coeff");
ADD_PROPERTY(m_C[3], "C3", FEProperty::Optional)->SetLongName("4th cv virial coeff");

END_FECORE_CLASS();

FERealGas::FERealGas(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_nva = 0;
    m_nvc = 1;
    m_R = m_Pr = m_Tr = 0;
    m_a0 = nullptr;
    m_A[0] = m_A[1] = m_A[2] = m_A[3] = m_A[4] = m_A[5] = m_A[6] = nullptr;
    m_C[0] = m_C[1] = m_C[2] = m_C[3] = nullptr;
}

//-----------------------------------------------------------------------------
//! initialization
bool FERealGas::Init()
{
    m_R  = GetGlobalConstant("R");
    m_Tr = GetGlobalConstant("T");
    m_Pr = GetGlobalConstant("P");
    
    if (m_R  <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section");    return false; }
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    if (m_Pr <= 0) { feLogError("A positive referential absolute pressure P must be defined in Globals section"); return false; }

    m_pMat = dynamic_cast<FEThermoFluid*>(GetParent());
    m_rhor = m_pMat->ReferentialDensity();

    // check if we should assume ideal gas
    if (m_nva == 0) {
		m_A[0] = fecore_alloc(FEConstFunction, GetFEModel());
		m_A[0]->SetParameter("value", 1.0);
        m_nva = 1;
    }
    m_a0->Init();
    for (int k=0; k<m_nva; ++k)
        if (m_A[k]) m_A[k]->Init();
    for (int k=0; k<m_nvc; ++k)
        if (m_C[k]) m_C[k]->Init();

    return true;
}

//-----------------------------------------------------------------------------
void FERealGas::Serialize(DumpStream& ar)
{
    FEElasticFluid::Serialize(ar);

    if (ar.IsShallow()) return;
    ar & m_pMat;
    ar & m_R & m_Pr & m_Tr & m_rhor;
}

//-----------------------------------------------------------------------------
//! gage pressure
double FERealGas::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double x = That/J;
    double y = x - 1;
    double p = 0;
    for (int k=1; k<=m_nva; ++k)
        p += m_A[k-1]->value(That)*pow(y,k);
    
    return p*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FERealGas::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double x = That/J;
    double y = x - 1;
    double dpJ = m_A[0]->value(That);
    for (int k=2; k<=m_nva; ++k)
        dpJ += k*m_A[k-1]->value(That)*pow(y,k-1);

    return -dpJ*m_Pr*That/pow(J,2);
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FERealGas::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double x = That/J;
    double y = x - 1;
    double dpJ2 = 2*m_A[0]->value(That);
    for (int k=2; k<=m_nva; ++k)
        dpJ2 += k*m_A[k-1]->value(That)*(2*y+(k-1)*That/J)*pow(y,k-2);

    return dpJ2*m_Pr*That/pow(J,3);
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FERealGas::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double y = That/J - 1;
    double dpT = 0;
    for (int k=1; k<=m_nva; ++k)
        dpT += (m_A[k-1]->derive(That)*y + k/J*m_A[k-1]->value(That))*pow(y,k-1);

    return dpT*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FERealGas::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double y = That/J - 1;
    double dpT2 = m_A[0]->deriv2(That)*y + 2/J*m_A[0]->derive(That);
    for (int k=2; k<=m_nva; ++k)
        dpT2 += (m_A[k-1]->deriv2(That)*pow(y,2) + 2*k/J*m_A[k-1]->derive(That)*y
                 + k*(k-1)/pow(J,2)*m_A[k-1]->value(That))*pow(y,k-2);

    return dpT2*m_Pr/pow(m_Tr,2);
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FERealGas::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double x = That/J;
    double y = x - 1;
    double dpJT = m_A[0]->value(That) + That*m_A[0]->derive(That);
    for (int k=2; k<=m_nva; ++k)
        dpJT += k*(That*m_A[k-1]->derive(That)*y + m_A[k-1]->value(That)*(k*x-1))*pow(y,k-2);

    return -dpJT*m_Pr/(m_Tr*pow(J, 2));
}

//-----------------------------------------------------------------------------
//! specific free energy
double FERealGas::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double x = That/J;
    double y = x - 1;
    double x2 = x*x, x3 = x2*x, x4 = x3*x, x5 = x4*x, x6 = x5*x;
    double lnx = log(x);
    double A[MAX_NVA], f[MAX_NVA];
    f[0] =  1 + x*(lnx-1);
    f[1] = -1 + x*(x-2*lnx);
    f[2] =  1 + x*(1.5+3*lnx) - 3*x2 + x3/2;
    f[3] = -1 - x*(10./3.+4*lnx) + 6*x2 - 2*x3 + x4/3;
    f[4] =  1 + x*(65./12.+5*lnx) - 10*x2 + 5*x3 - 5*x4/3 + x5/4;
    f[5] = -1 - x*(77./10.+6*lnx) + 15*x2 - 10*x3 + 5*x4 - 1.5*x5 + x6/5;
    f[6] =  1 + x*(203./20.+7*lnx) - 21*x2 + 35*x3/2 - 35*x4/3 + 21*x5/4 - 7*x6/5;
    for (int k=0; k<m_nva; ++k)
        A[k] = m_A[k]->value(That);
    double a = m_a0->value(That);
    for (int k=0; k<m_nva; ++k)
        a += A[k]*J*f[k];

    return a*m_Pr/m_rhor;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FERealGas::SpecificEntropy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double x = That/J;
    double y = x - 1;
    double x2 = x*x, x3 = x2*x, x4 = x3*x, x5 = x4*x, x6 = x5*x;
    double lnx = log(x);
    double A[MAX_NVA], dA[MAX_NVA], f[MAX_NVA], df[MAX_NVA];
    f[0] =  1 + x*(lnx-1);
    f[1] = -1 + x*(x-2*lnx);
    f[2] =  1 + x*(1.5+3*lnx) - 3*x2 + x3/2;
    f[3] = -1 - x*(10./3.+4*lnx) + 6*x2 - 2*x3 + x4/3;
    f[4] =  1 + x*(65./12.+5*lnx) - 10*x2 + 5*x3 - 5*x4/3 + x5/4;
    f[5] = -1 - x*(77./10.+6*lnx) + 15*x2 - 10*x3 + 5*x4 - 1.5*x5 + x6/5;
    f[6] =  1 + x*(203./20.+7*lnx) - 21*x2 + 35*x3/2 - 35*x4/3 + 21*x5/4 - 7*x6/5;
    df[0] = lnx;
    df[1] = 2*(-1 + x - lnx);
    df[2] = 1.5*(3 - 4*x + x2 + 2*lnx);
    df[3] = 2./3.*(-11 + 18*x - 9*x2 + 2*x3 - 6*lnx);
    df[4] = 5./12.*(25 - 48*x + 36*x2 - 16*x3 + 3*x4 + 12*lnx);
    df[5] = -13.7 + 30*x - 30*x2 + 20*x3 - 7.5*x4 + 1.2*x5 - 6*lnx;
    df[6] = 7./60.*(147 - 360*x + 450*x2 - 400*x3 + 225*x4 - 72*x5 + 10*x6 + 60*lnx);
    for (int k=0; k<m_nva; ++k) {
        A[k] = m_A[k]->value(That);
        dA[k] = m_A[k]->derive(That);
    }
    double s = -m_a0->derive(That);
    for (int k=0; k<m_nva; ++k)
        s -= A[k]*df[k] + J*dA[k]*f[k];

    return s*m_Pr/(m_Tr*m_rhor);
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FERealGas::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    // get the specific free energy
    double a = SpecificFreeEnergy(mp);
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double a0 = m_a0->value(That)*m_Pr/m_rhor;
    
    // the specific strain energy is the difference between these two values
    return a - a0;
}

//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FERealGas::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();

    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double y = That/J - 1;
    double cv = m_C[0]->value(That);
    for (int k=1; k<m_nvc; ++k)
        cv += m_C[k]->value(That)*pow(y,k);

    return cv*m_Pr/(m_Tr*m_rhor);
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FERealGas::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double x = That/J;
    double y = x - 1;
    double dcvJ = 0;
    for (int k=1; k<m_nvc; ++k)
        dcvJ -= m_C[k]->value(That)*pow(y,k-1);
    return dcvJ*m_Pr/(m_Tr*m_rhor)*x/J;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FERealGas::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double x = That/J;
    double y = x - 1;
    double dcvT = m_C[0]->derive(That);
    for (int k=1; k<m_nvc; ++k)
        dcvT += (m_C[k]->derive(That)*y+k/J*m_C[k]->value(That))*pow(y,k-1);
    return dcvT*m_Pr/(pow(m_Tr,2)*m_rhor);
}

//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FERealGas::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
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
bool FERealGas::Dilatation(const double T, const double p, double& e)
{
    double errrel = 1e-6;
    double errabs = 1e-6;
    int maxiter = 100;
    switch (m_nva) {
        case 0:
            return false;
            break;
        case 1:
        {
            double q = T/m_Tr;
            double x = p/m_Pr/m_A[0]->value(q);
            e = (q+1)/(x+1) - 1;
            return true;
        }
            break;
        case 2:
        {
            double q = T/m_Tr;
            double A1 = m_A[0]->value(q);
            double A2 = m_A[1]->value(q);
            double det = A1*A1 + 4*A2*p/m_Pr;
            if (det < 0) return false;
            det = sqrt(det);
            // only one root of this quadratic equation is valid.
            double x = (-A1+det)/(2*A2);
            e = (q+1)/(x+1) - 1;
            return true;
        }
            break;
            
        default:
        {
            FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
            FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
            ft->m_T = T;
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
    return false;
}
