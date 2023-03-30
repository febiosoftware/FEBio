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
ADD_PARAMETER(m_nvc  , FE_RANGE_GREATER_OR_EQUAL(0), "nvc");

ADD_PROPERTY(m_a0  , "a0");
ADD_PROPERTY(m_A[0], "A1", FEProperty::Optional);
ADD_PROPERTY(m_A[1], "A2", FEProperty::Optional);
ADD_PROPERTY(m_A[2], "A3", FEProperty::Optional);
ADD_PROPERTY(m_A[3], "A4", FEProperty::Optional);
ADD_PROPERTY(m_A[4], "A5", FEProperty::Optional);
ADD_PROPERTY(m_A[5], "A6", FEProperty::Optional);
ADD_PROPERTY(m_A[6], "A7", FEProperty::Optional);

END_FECORE_CLASS();

FERealGas::FERealGas(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_nvc = 0;
    m_R = m_Pr = m_Tr = 0;
    m_a0 = nullptr;
    m_A[0] = m_A[1] = m_A[2] = m_A[3] = m_A[4] = m_A[5] = m_A[6] = nullptr;
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
    if (m_nvc == 0) {
		m_A[0] = fecore_alloc(FEConstFunction, GetFEModel());
		m_A[0]->SetParameter("value", 1.0);
        m_nvc = 1;
    }
    m_a0->Init();
    for (int k=0; k<m_nvc; ++k)
        if (m_A[k]) m_A[k]->Init();
    
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
    
    double dq = tf.m_T/m_Tr;
    double J = 1 + fp.m_ef;
    double x = (1+dq)/J - 1;
    double p = 0;
    for (int k=1; k<=m_nvc; ++k)
        p += m_A[k-1]->value(dq)*pow(x,k);
    
    return p*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FERealGas::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double dq = tf.m_T/m_Tr;
    double q = 1 + dq;
    double J = 1 + fp.m_ef;
    double x = q/J - 1;
    double dpJ = m_A[0]->value(dq);
    for (int k=2; k<=m_nvc; ++k)
        dpJ += k*m_A[k-1]->value(dq)*pow(x,k-1);

    return -dpJ*m_Pr*q/pow(J,2);
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FERealGas::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double dq = tf.m_T/m_Tr;
    double q = 1 + dq;
    double J = 1 + fp.m_ef;
    double x = q/J - 1;
    double dpJ2 = 2*m_A[0]->value(dq);
    for (int k=2; k<=m_nvc; ++k)
        dpJ2 += k*m_A[k-1]->value(dq)*(2*x+(k-1)*q/J)*pow(x,k-2);

    return dpJ2*m_Pr*q/pow(J,3);
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FERealGas::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double dq = tf.m_T/m_Tr;
    double q = 1 + dq;
    double J = 1 + fp.m_ef;
    double x = q/J - 1;
    double dpT = 0;
    for (int k=1; k<=m_nvc; ++k)
        dpT += (m_A[k-1]->derive(dq)*x + k/J*m_A[k-1]->value(dq))*pow(x,k-1);

    return dpT*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FERealGas::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double dq = tf.m_T/m_Tr;
    double q = 1 + dq;
    double J = 1 + fp.m_ef;
    double x = q/J - 1;
    double dpT2 = m_A[0]->deriv2(dq)*x + 2/J*m_A[0]->derive(dq);
    for (int k=2; k<=m_nvc; ++k)
        dpT2 += (m_A[k-1]->deriv2(dq)*pow(x,2) + 2*k/J*m_A[k-1]->derive(dq)*x
                 + k*(k-1)/pow(J,2)*m_A[k-1]->value(dq))*pow(x,k-2);

    return dpT2*m_Pr/pow(m_Tr,2);
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FERealGas::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double dq = tf.m_T/m_Tr;
    double q = 1 + dq;
    double J = 1 + fp.m_ef;
    double x = q/J - 1;
    double dpJT = m_A[0]->value(dq) + q*m_A[0]->derive(dq);
    for (int k=2; k<=m_nvc; ++k)
        dpJT += k*(q*m_A[k-1]->derive(dq)*x + m_A[k-1]->value(dq)*(x+(k-1)*q/J))*pow(x,k-2);

    return -dpJT*m_Pr/(m_Tr*pow(J, 2));
}

//-----------------------------------------------------------------------------
//! specific free energy
double FERealGas::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double dq = tf.m_T/m_Tr;
    double q = 1 + dq;
    double q2 = q*q, q3 = q2*q, q4 = q3*q, q5 = q4*q, q6 = q5*q, q7 = q6*q;
    double lnq = log(q);
    double J = 1 + fp.m_ef;
    double J2 = J*J, J3 = J2*J, J4 = J3*J, J5 = J4*J, J6 = J5*J, J7 = J6*J;
    double lnJ = log(J);
    double A[MAX_NVC];
    for (int k=0; k<m_nvc; ++k)
        A[k] = m_A[k]->value(dq);
    double a = m_a0->value(dq);
    a += (-60*J6*q*A[0] - 60*J6*lnJ*q*A[0] + 60*J6*lnq*q*A[0] + 120*J6*lnJ*q*A[1] -
          120*J6*lnq*q*A[1] + 60*J5*q2*A[1] + 90*J6*q*A[2] - 180*J6*lnJ*q*A[2] +
          180*J6*lnq*q*A[2] - 180*J5*q2*A[2] + 30*J4*q3*A[2] - 200*J6*q*A[3] +
          240*J6*lnJ*q*A[3] - 240*J6*lnq*q*A[3] + 360*J5*q2*A[3] - 120*J4*q3*A[3] +
          20*J3*q4*A[3] + 325*J6*q*A[4] - 300*J6*lnJ*q*A[4] + 300*J6*lnq*q*A[4] -
          600*J5*q2*A[4] + 300*J4*q3*A[4] - 100*J3*q4*A[4] + 15*J2*q5*A[4] -
          462*J6*q*A[5] + 360*J6*lnJ*q*A[5] - 360*J6*lnq*q*A[5] + 900*J5*q2*A[5] -
          600*J4*q3*A[5] + 300*J3*q4*A[5] - 90*J2*q5*A[5] + 12*J*q6*A[5] + 609*J6*q*A[6] -
          420*J6*lnJ*q*A[6] + 420*J6*lnq*q*A[6] - 1260*J5*q2*A[6] + 1050*J4*q3*A[6] -
          700*J3*q4*A[6] + 315*J2*q5*A[6] - 84*J*q6*A[6] + 10*q7*A[6] +
          60*J7*(A[0] - A[1] + A[2] - A[3] + A[4] - A[5] + A[6]))/
    (60.*J6);

    return a*m_Pr/m_rhor;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FERealGas::SpecificEntropy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double dq = tf.m_T/m_Tr;
    double q = 1 + dq;
    double q2 = q*q, q3 = q2*q, q4 = q3*q, q5 = q4*q, q6 = q5*q, q7 = q6*q;
    double lnq = log(q);
    double J = 1 + fp.m_ef;
    double J2 = J*J, J3 = J2*J, J4 = J3*J, J5 = J4*J, J6 = J5*J, J7 = J6*J;
    double lnJ = log(J);
    double A[MAX_NVC], dA[MAX_NVC];
    for (int k=0; k<m_nvc; ++k) {
        A[k] = m_A[k]->value(dq);
        dA[k] = m_A[k]->derive(dq);
    }
    double da0 = m_a0->derive(dq);
    double s = -((60*da0*J6 + 120*J5*q*A[1] - 360*J5*q*A[2] + 90*J4*q2*A[2] + 720*J5*q*A[3] -
                       360*J4*q2*A[3] + 80*J3*q3*A[3] - 1200*J5*q*A[4] + 900*J4*q2*A[4] -
                       400*J3*q3*A[4] + 75*J2*q4*A[4] + 1800*J5*q*A[5] - 1800*J4*q2*A[5] +
                       1200*J3*q3*A[5] - 450*J2*q4*A[5] + 72*J*q5*A[5] - 2520*J5*q*A[6] +
                       3150*J4*q2*A[6] - 2800*J3*q3*A[6] + 1575*J2*q4*A[6] - 504*J*q5*A[6] +
                       70*q6*A[6] + 60*J7*dA[0] - 60*J7*dA[1] + 60*J5*q2*dA[1] + 60*J7*dA[2] -
                       180*J5*q2*dA[2] + 30*J4*q3*dA[2] - 60*J7*dA[3] + 360*J5*q2*dA[3] -
                       120*J4*q3*dA[3] + 20*J3*q4*dA[3] + 60*J7*dA[4] - 600*J5*q2*dA[4] +
                       300*J4*q3*dA[4] - 100*J3*q4*dA[4] + 15*J2*q5*dA[4] - 60*J7*dA[5] +
                       900*J5*q2*dA[5] - 600*J4*q3*dA[5] + 300*J3*q4*dA[5] - 90*J2*q5*dA[5] +
                       12*J*q6*dA[5] + 60*J7*dA[6] - 1260*J5*q2*dA[6] + 1050*J4*q3*dA[6] -
                       700*J3*q4*dA[6] + 315*J2*q5*dA[6] - 84*J*q6*dA[6] + 10*q7*dA[6] +
                       J6*(-120*A[1] + 270*A[2] - 440*A[3] + 625*A[4] - 822*A[5] + 1029*A[6] -
                           60*q*dA[0] + 90*q*dA[2] - 200*q*dA[3] + 325*q*dA[4] - 462*q*dA[5] +
                           609*q*dA[6] - 60*lnJ*(A[0] - 2*A[1] + 3*A[2] - 4*A[3] + 5*A[4] - 6*A[5] +
                                                 7*A[6] + q*dA[0] - 2*q*dA[1] + 3*q*dA[2] - 4*q*dA[3] + 5*q*dA[4] -
                                                 6*q*dA[5] + 7*q*dA[6]) +
                           60*lnq*(A[0] - 2*A[1] + 3*A[2] - 4*A[3] + 5*A[4] - 6*A[5] + 7*A[6] +
                                   q*dA[0] - 2*q*dA[1] + 3*q*dA[2] - 4*q*dA[3] + 5*q*dA[4] - 6*q*dA[5] +
                                   7*q*dA[6]))))/(60.*J6);

    return s*m_Pr/(m_Tr*m_rhor);
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FERealGas::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    // get the specific free energy
    double a = SpecificFreeEnergy(mp);
    double dq = tf.m_T/m_Tr;
    double a0 = m_a0->value(dq)*m_Pr/m_rhor;
    
    // the specific strain energy is the difference between these two values
    return a - a0;
}

//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FERealGas::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();

    double dq = tf.m_T/m_Tr;
    double q = 1 + dq;
    double q2 = q*q, q3 = q2*q, q4 = q3*q, q5 = q4*q, q6 = q5*q, q7 = q6*q;
    double lnq = log(q);
    double J = 1 + fp.m_ef;
    double J2 = J*J, J3 = J2*J, J4 = J3*J, J5 = J4*J, J6 = J5*J, J7 = J6*J;
    double lnJ = log(J);
    double A[MAX_NVC], dA[MAX_NVC], d2A[MAX_NVC];
    for (int k=0; k<m_nvc; ++k) {
        A[k] = m_A[k]->value(dq);
        dA[k] = m_A[k]->derive(dq);
        d2A[k] = m_A[k]->deriv2(dq);
    }
    double d2a0 = m_a0->deriv2(dq);
    double cv = -q*(60*d2a0*J6 + 120*J5*A[1] - (120*J*J5*A[1])/q - 360*J5*A[2] + 180*J4*q*A[2] +
                    720*J5*A[3] - 720*J4*q*A[3] + 240*J3*q2*A[3] - 1200*J5*A[4] + 1800*J4*q*A[4] -
                    1200*J3*q2*A[4] + 300*J2*q3*A[4] + 1800*J5*A[5] - 3600*J4*q*A[5] +
                    3600*J3*q2*A[5] - 1800*J2*q3*A[5] + 360*J*q4*A[5] - 2520*J5*A[6] +
                    6300*J4*q*A[6] - 8400*J3*q2*A[6] + 6300*J2*q3*A[6] - 2520*J*q4*A[6] +
                    420*q5*A[6] + 60*J7*d2A[0] - 60*J7*d2A[1] + 60*J5*q2*d2A[1] + 60*J7*d2A[2] -
                    180*J5*q2*d2A[2] + 30*J4*q3*d2A[2] - 60*J7*d2A[3] + 360*J5*q2*d2A[3] -
                    120*J4*q3*d2A[3] + 20*J3*q4*d2A[3] + 60*J7*d2A[4] - 600*J5*q2*d2A[4] +
                    300*J4*q3*d2A[4] - 100*J3*q4*d2A[4] + 15*J2*q5*d2A[4] - 60*J7*d2A[5] +
                    900*J5*q2*d2A[5] - 600*J4*q3*d2A[5] + 300*J3*q4*d2A[5] - 90*J2*q5*d2A[5] +
                    12*J*q6*d2A[5] + 60*J7*d2A[6] - 1260*J5*q2*d2A[6] + 1050*J4*q3*d2A[6] -
                    700*J3*q4*d2A[6] + 315*J2*q5*d2A[6] - 84*J*q6*d2A[6] + 10*q7*d2A[6] +
                    240*J5*q*dA[1] - 720*J5*q*dA[2] + 180*J4*q2*dA[2] + 1440*J5*q*dA[3] -
                    720*J4*q2*dA[3] + 160*J3*q3*dA[3] - 2400*J5*q*dA[4] + 1800*J4*q2*dA[4] -
                    800*J3*q3*dA[4] + 150*J2*q4*dA[4] + 3600*J5*q*dA[5] - 3600*J4*q2*dA[5] +
                    2400*J3*q3*dA[5] - 900*J2*q4*dA[5] + 144*J*q5*dA[5] - 5040*J5*q*dA[6] +
                    6300*J4*q2*dA[6] - 5600*J3*q3*dA[6] + 3150*J2*q4*dA[6] - 1008*J*q5*dA[6] +
                    140*q6*dA[6] + J6*((60*(A[0] + 3*A[2] - 4*A[3] + 5*A[4] - 6*A[5] + 7*A[6]))/q +
                                       q*(-60*(1 + lnJ - lnq)*d2A[0] - 120*lnq*d2A[1] + 90*d2A[2] +
                                          180*lnq*d2A[2] - 200*d2A[3] - 240*lnq*d2A[3] + 325*d2A[4] +
                                          300*lnq*d2A[4] - 462*d2A[5] - 360*lnq*d2A[5] +
                                          60*lnJ*(2*d2A[1] - 3*d2A[2] + 4*d2A[3] - 5*d2A[4] + 6*d2A[5] -
                                                  7*d2A[6]) + 609*d2A[6] + 420*lnq*d2A[6]) -
                                       2*(120*dA[1] - 270*dA[2] + 440*dA[3] - 625*dA[4] + 822*dA[5] - 1029*dA[6] +
                                          60*lnJ*(dA[0] - 2*dA[1] + 3*dA[2] - 4*dA[3] + 5*dA[4] - 6*dA[5] +
                                                  7*dA[6]) - 60*lnq*(dA[0] - 2*dA[1] + 3*dA[2] - 4*dA[3] + 5*dA[4] -
                                                                     6*dA[5] + 7*dA[6]))))/(60.*J6);

    return cv*m_Pr/(m_Tr*m_rhor);
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FERealGas::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double dq = tf.m_T/m_Tr;
    double q = 1 + dq;
    double q2 = q*q, q3 = q2*q, q4 = q3*q, q5 = q4*q, q6 = q5*q, q7 = q6*q;
    double lnq = log(q);
    double J = 1 + fp.m_ef;
    double J2 = J*J, J3 = J2*J, J4 = J3*J, J5 = J4*J, J6 = J5*J, J7 = J6*J;
    double lnJ = log(J);
    double A[MAX_NVC], dA[MAX_NVC], d2A[MAX_NVC];
    for (int k=0; k<m_nvc; ++k) {
        A[k] = m_A[k]->value(dq);
        dA[k] = m_A[k]->derive(dq);
        d2A[k] = m_A[k]->deriv2(dq);
    }
    double d2a0 = m_a0->deriv2(dq);
    double dcvJ = -(q*(60*d2a0*J6 + 120*J5*A[1] - (120*J*J5*A[1])/q - 360*J5*A[2] + 180*J4*q*A[2] +
                            720*J5*A[3] - 720*J4*q*A[3] + 240*J3*q2*A[3] - 1200*J5*A[4] + 1800*J4*q*A[4] -
                            1200*J3*q2*A[4] + 300*J2*q3*A[4] + 1800*J5*A[5] - 3600*J4*q*A[5] +
                            3600*J3*q2*A[5] - 1800*J2*q3*A[5] + 360*J*q4*A[5] - 2520*J5*A[6] +
                            6300*J4*q*A[6] - 8400*J3*q2*A[6] + 6300*J2*q3*A[6] - 2520*J*q4*A[6] +
                            420*q5*A[6] + 60*J7*d2A[0] - 60*J7*d2A[1] + 60*J5*q2*d2A[1] + 60*J7*d2A[2] -
                            180*J5*q2*d2A[2] + 30*J4*q3*d2A[2] - 60*J7*d2A[3] + 360*J5*q2*d2A[3] -
                            120*J4*q3*d2A[3] + 20*J3*q4*d2A[3] + 60*J7*d2A[4] - 600*J5*q2*d2A[4] +
                            300*J4*q3*d2A[4] - 100*J3*q4*d2A[4] + 15*J2*q5*d2A[4] - 60*J7*d2A[5] +
                            900*J5*q2*d2A[5] - 600*J4*q3*d2A[5] + 300*J3*q4*d2A[5] - 90*J2*q5*d2A[5] +
                            12*J*q6*d2A[5] + 60*J7*d2A[6] - 1260*J5*q2*d2A[6] + 1050*J4*q3*d2A[6] -
                            700*J3*q4*d2A[6] + 315*J2*q5*d2A[6] - 84*J*q6*d2A[6] + 10*q7*d2A[6] +
                            240*J5*q*dA[1] - 720*J5*q*dA[2] + 180*J4*q2*dA[2] + 1440*J5*q*dA[3] -
                            720*J4*q2*dA[3] + 160*J3*q3*dA[3] - 2400*J5*q*dA[4] + 1800*J4*q2*dA[4] -
                            800*J3*q3*dA[4] + 150*J2*q4*dA[4] + 3600*J5*q*dA[5] - 3600*J4*q2*dA[5] +
                            2400*J3*q3*dA[5] - 900*J2*q4*dA[5] + 144*J*q5*dA[5] - 5040*J5*q*dA[6] +
                            6300*J4*q2*dA[6] - 5600*J3*q3*dA[6] + 3150*J2*q4*dA[6] - 1008*J*q5*dA[6] +
                            140*q6*dA[6] + J6*((60*(A[0] + 3*A[2] - 4*A[3] + 5*A[4] - 6*A[5] + 7*A[6]))/q +
                                               q*(-60*(1 + lnJ - lnq)*d2A[0] - 120*lnq*d2A[1] + 90*d2A[2] +
                                                  180*lnq*d2A[2] - 200*d2A[3] - 240*lnq*d2A[3] + 325*d2A[4] +
                                                  300*lnq*d2A[4] - 462*d2A[5] - 360*lnq*d2A[5] +
                                                  60*lnJ*(2*d2A[1] - 3*d2A[2] + 4*d2A[3] - 5*d2A[4] + 6*d2A[5] -
                                                          7*d2A[6]) + 609*d2A[6] + 420*lnq*d2A[6]) -
                                               2*(120*dA[1] - 270*dA[2] + 440*dA[3] - 625*dA[4] + 822*dA[5] - 1029*dA[6] +
                                                  60*lnJ*(dA[0] - 2*dA[1] + 3*dA[2] - 4*dA[3] + 5*dA[4] - 6*dA[5] +
                                                          7*dA[6]) - 60*lnq*(dA[0] - 2*dA[1] + 3*dA[2] - 4*dA[3] + 5*dA[4] -
                                                                             6*dA[5] + 7*dA[6])))))/(60.*J6);
    
    return dcvJ*m_Pr/(m_Tr*m_rhor);
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FERealGas::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = m_Tr + tf.m_T;
    double dcvT = IsochoricSpecificHeatCapacity(mp)/T;  // this is incomplete
    return dcvT;
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
bool FERealGas::Dilatation(const double T, const double p, const double c, double& e)
{
    double errrel = 1e-6;
    double errabs = 1e-6;
    int maxiter = 100;
    switch (m_nvc) {
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
