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
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include "FEFluidMaterialPoint.h"
#include "FEThermoFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERealLiquid, FEElasticFluid)

// material parameters
ADD_PARAMETER(m_nvc  , FE_RANGE_GREATER(0), "nvc");

ADD_PROPERTY(m_psat, "psat");
ADD_PROPERTY(m_asat, "asat");
ADD_PROPERTY(m_ssat, "ssat");
ADD_PROPERTY(m_Jsat, "Jsat");
ADD_PROPERTY(m_B[0], "B1"  );
ADD_PROPERTY(m_B[1], "B2", FEProperty::Optional);
ADD_PROPERTY(m_B[2], "B3", FEProperty::Optional);

END_FECORE_CLASS();

FERealLiquid::FERealLiquid(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_R = m_Pr = m_Tr = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FERealLiquid::Init()
{
    m_R  = GetFEModel()->GetGlobalConstant("R");
    m_Tr = GetFEModel()->GetGlobalConstant("T");
    m_Pr = GetFEModel()->GetGlobalConstant("P");
    
    if (m_R  <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section");    return false; }
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    if (m_Pr <= 0) { feLogError("A positive referential absolute pressure P must be defined in Globals section"); return false; }
    
    m_pMat = dynamic_cast<FEThermoFluid*>(GetParent());
    m_rhor = m_pMat->ReferentialDensity();
    
    return true;
}

//-----------------------------------------------------------------------------
//! gage pressure
double FERealLiquid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double q = tf.m_T/m_Tr;
    double p = m_psat->value(q);
    double de = fp.m_Jf - m_Jsat->value(q);
    for (int k=1; k<=m_nvc; ++k)
        p += m_B[k-1]->value(q)*pow(de,k);
    
    return p*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FERealLiquid::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double q = tf.m_T/m_Tr;
    double dpJ = m_B[0]->value(q);
    double de = fp.m_Jf - m_Jsat->value(q);
    for (int k=2; k<=m_nvc; ++k)
        dpJ += k*m_B[k-1]->value(q)*pow(de,k-1);
    
    return dpJ*m_Pr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FERealLiquid::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double q = tf.m_T/m_Tr;
    double dpJ2 = 2*m_B[1]->value(q);
    double de = fp.m_Jf - m_Jsat->value(q);
    for (int k=3; k<=m_nvc; ++k)
        dpJ2 += k*(k-1)*m_B[k-1]->value(q)*pow(de,k-2);
    
    return dpJ2*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FERealLiquid::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double q = tf.m_T/m_Tr;
    double dpT = m_psat->derive(q);
    double dJsat = m_Jsat->derive(q);
    double de = fp.m_Jf - m_Jsat->value(q);
    for (int k=1; k<=m_nvc; ++k)
        dpT += m_B[k-1]->derive(q)*pow(de,k)
        - k*dJsat*m_B[k-1]->value(q)*pow(de,k-1);
    
    return dpT*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FERealLiquid::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double q = tf.m_T/m_Tr;
    double dpT2 = m_psat->deriv2(q);
    double dJsat = m_Jsat->derive(q);
    double dJsa2 = m_Jsat->deriv2(q);
    double de = fp.m_Jf - m_Jsat->value(q);
    vector<double> B(m_nvc,0);
    for (int k=0; k<m_nvc; ++k) B[k] = m_B[k]->value(q);
    for (int k=1; k<=m_nvc; ++k)
        dpT2 += m_B[k-1]->deriv2(q)*pow(de,k)
        - (2*dJsat*m_B[k-1]->derive(q) + dJsa2*B[k-1])*k*pow(de,k-1);
    for (int k=2; k<=m_nvc; ++k)
        dpT2 += k*(k-1)*pow(dJsat,2)*B[k-1]*pow(de,k-2);

    return dpT2*m_Pr/pow(m_Tr,2);
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FERealLiquid::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double q = tf.m_T/m_Tr;
    double dpJT = m_B[0]->derive(q);
    double dJsat = m_Jsat->derive(q);
    double de = fp.m_Jf - m_Jsat->value(q);
    for (int k=2; k<=m_nvc; ++k)
        dpJT += k*m_B[k-1]->derive(q)*pow(de,k-1)
        - k*(k-1)*dJsat*m_B[k-1]->value(q)*pow(de,k-2);
    
    return dpJT*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FERealLiquid::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double q = tf.m_T/m_Tr;
    double de = fp.m_Jf - m_Jsat->value(q);
    double a = m_asat->value(q) - de*m_psat->value(q);
    for (int k=1; k<=m_nvc; ++k)
        a -= m_B[k-1]->value(q)/(k+1)*pow(de,k+1);
    
    return a*m_Pr/m_rhor;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FERealLiquid::SpecificEntropy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double q = tf.m_T/m_Tr;
    double de = fp.m_Jf - m_Jsat->value(q);
    double dJsat = m_Jsat->derive(q);
    double s = m_ssat->value(q) - m_psat->value(q)*dJsat + de*m_psat->derive(q);
    for (int k=1; k<=m_nvc; ++k)
        s += (m_B[k-1]->derive(q)*de/(k+1) - dJsat*m_B[k-1]->value(q))*pow(de,k);

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
    
    // evaluate the dilatation that makes the pressure = 0 at the given temperature
    double e = fp.m_Jf - 1;
    bool good = Dilatation(tf.m_T, 0, e);
    assert(good);
    
    // for this dilatation evaluate the specific free energy
    FEFluidMaterialPoint* fmp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* fmt = new FEThermoFluidMaterialPoint(fmp);
    fmp->m_Jf = 1 + e;
    fmt->m_T = tf.m_T;
    double a0 = SpecificFreeEnergy(*fmt);

    // the specific strain energy is the difference between these two values
    return a - a0;
}

//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FERealLiquid::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double q = tf.m_T/m_Tr;
    double de = fp.m_Jf - m_Jsat->value(q);
    double dJsat = m_Jsat->derive(q);
    double dJsa2 = m_Jsat->deriv2(q);
    double dsT = m_ssat->derive(q) + de*m_psat->deriv2(q)
    - 2*dJsat*m_psat->derive(q) - dJsa2*m_psat->value(q);
    for (int k=1; k<=m_nvc; ++k)
        dsT += (k*pow(dJsat, 2) - de*dJsa2)*m_B[k-1]->value(q)*pow(de, k-1)
        - 2*dJsat*m_B[k-1]->derive(q)*pow(de, k)
        + m_B[k-1]->deriv2(q)/(k+1)*pow(de, k+1);
    
    return (m_Tr + tf.m_T)*dsT*m_Pr/(m_rhor*pow(m_Tr,2));
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FERealLiquid::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double q = tf.m_T/m_Tr;
    double de = fp.m_Jf - m_Jsat->value(q);
    double dJsat = m_Jsat->derive(q);
    double dJsa2 = m_Jsat->deriv2(q);
    double dsTJ = m_psat->deriv2(q) - dJsa2*m_B[0]->value(q);
    for (int k=1; k<=m_nvc; ++k)
        dsTJ += (de*m_B[k-1]->deriv2(q) - 2*dJsat*k*m_B[k-1]->derive(q))*pow(de, k-1);
    for (int k=2; k<=m_nvc; ++k)
        dsTJ += k*m_B[k-1]->value(q)*((k-1)*pow(dJsat, 2) - de*dJsa2)*pow(de, k-2);
    
    return (m_Tr + tf.m_T)*dsTJ*m_Pr/(m_rhor*pow(m_Tr,2));
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FERealLiquid::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = m_Tr + tf.m_T;
    double dcvT = IsochoricSpecificHeatCapacity(mp)/T;  // this is incomplete
    return dcvT;
}

//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FERealLiquid::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    // evaluate the dilatation that makes the gage pressure = 0 at the given temperature
    double e = fp.m_Jf - 1; // initial guess
    bool good = Dilatation(tf.m_T, 0, e);
    assert(good);

    // create a material point with this temperature and dilatation
    FEFluidMaterialPoint* fmp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* fmt = new FEThermoFluidMaterialPoint(fmp);
    fmp->m_Jf = 1 + e;
    fmt->m_T = tf.m_T;
    
    // for this dilatation evaluate the following terms
    double cv = IsochoricSpecificHeatCapacity(*fmt);
    double dpT = Tangent_Temperature(*fmt);
    double dpJ = Tangent_Strain(*fmt);
    double cp = cv - (m_Tr + tf.m_T)*pow(dpT, 2)/(m_rhor*dpJ);
    
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
    ft->m_T = T;
    double q = T/m_Tr;
    double psat = m_psat->value(q);
    double Jsat = m_Jsat->value(q);
    if (p == psat) { e = Jsat - 1; return true; }
    switch (m_nvc) {
        case 0:
            return false;
            break;
        case 1:
        {
            double B0 = m_B[0]->value(q);
            if (B0 == 0) { e = Jsat - 1; return true; }
            e = (p/m_Pr - psat)/B0 + Jsat - 1;
            return true;
        }
            break;
        case 2:
        {
            double B0 = m_B[0]->value(q);
            double B1 = m_B[1]->value(q);
            if ((B0 == 0) || (B1 == 0)) { e = Jsat - 1; return true; }
            double det = pow(B0,2) + 4*B1*(p/m_Pr - psat);
            if (det < 0) return false;
            det = sqrt(det);
            double e0 = Jsat - 1 - (B0+det)/(2*B1);
//            double e1 = Jsat - 1 - (B0-det)/(2*B1);
            e = e0;
            return true;
        }
            break;

        default:
        {
            bool convgd = false;
            bool done = false;
            int iter = 0;
            do {
                ++iter;
                fp->m_Jf = 1 + e;
                double f = Pressure(*ft) - p;
                double df = Tangent_Strain(*ft);
                double de = -f/df;
                e += de;
                if ((fabs(de) < errrel*fabs(e)) ||
                    (fabs(f/m_Pr) < errabs)) { convgd = true; done = true; }
                if (iter > maxiter) done = true;
            } while (!done);
            
            return convgd;
        }
            break;
    }
    return false;
}

//-----------------------------------------------------------------------------
//! pressure from state variables
double FERealLiquid::Pressure(const double ef, const double T)
{
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_Jf = 1 + ef;
    ft->m_T = T;
    return Pressure(*ft);
}
