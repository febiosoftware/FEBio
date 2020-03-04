/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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


#include "FEIdealGas.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include "FEFluidMaterialPoint.h"
#include "FEThermoFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEIdealGas, FEElasticFluid)

    // material parameters
    ADD_PARAMETER(m_M   , FE_RANGE_GREATER(0.0), "M");
    ADD_PARAMETER(m_ar  , "ar");
    ADD_PARAMETER(m_sr  , "sr");
    ADD_PARAMETER(m_C[0], "C0");
    ADD_PARAMETER(m_C[1], "C1");
    ADD_PARAMETER(m_C[2], "C2");
    ADD_PARAMETER(m_C[3], "C3");
    ADD_PARAMETER(m_C[4], "C4");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! initialization
bool FEIdealGas::Init()
{
    m_R  = GetFEModel()->GetGlobalConstant("R");
    m_Tr = GetFEModel()->GetGlobalConstant("T");
    m_Pr = GetFEModel()->GetGlobalConstant("P");
    
    if (m_R  <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section");    return false; }
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    if (m_Pr <= 0) { feLogError("A positive referential absolute pressure P must be defined in Globals section"); return false; }
    
    return true;
}

//-----------------------------------------------------------------------------
//! gage pressure
double FEIdealGas::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = m_Tr + tf.m_T;
    
    double p = m_Pr*(T/(fp.m_Jf*m_Tr) - 1);

    return p;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FEIdealGas::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = m_Tr + tf.m_T;

    double dp = -m_Pr*T/(fp.m_Jf*fp.m_Jf*m_Tr);

    return dp;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FEIdealGas::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = m_Tr + tf.m_T;

    double d2p = 2*m_Pr*T/(pow(fp.m_Jf,3)*m_Tr);

    return d2p;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FEIdealGas::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    
    double dp = m_Pr/(fp.m_Jf*m_Tr);

    return dp;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FEIdealGas::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FEIdealGas::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    
    double d2p = -m_Pr/(fp.m_Jf*fp.m_Jf*m_Tr);

    return d2p;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FEIdealGas::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = fp.m_Jf;
    double T = tf.m_T + m_Tr;
    double Tr = m_Tr;
    
    // referential free energy
    double a = m_ar - m_sr*(T-Tr);
    
    // add a circle
    a += m_R/m_M*(m_C[0]*(T-Tr-T*log(T/Tr))
                  +pow(T-Tr, 2)/2*(m_C[1]+m_C[2]/3*(T+2*Tr)
                                   +m_C[3]/6*(T*T+2*Tr*T+3*Tr*Tr)
                                   +m_C[4]/10*(pow(T, 3)+2*Tr*T*T+3*Tr*Tr*T+4*pow(Tr, 3))));
    
    // add strain-dependent contribution
    a += m_R/m_M*(J*Tr-T+T*log(T/(J*Tr)));

    return a;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FEIdealGas::SpecificEntropy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = fp.m_Jf;
    double T = tf.m_T + m_Tr;
    double Tr = m_Tr;
    
    // referential entropy
    double s = m_sr;
    
    // add s circle
    s += m_R/m_M*(m_C[0]*log(T/Tr) + (T-Tr)*(m_C[1]+m_C[3]/3*(T*T+Tr*T+Tr*Tr))
                  +(T*T-Tr*Tr)/2*(m_C[2]+m_C[4]/2*(T*T+Tr*Tr)));
    
    // add strain-dependent contribution
    s += -m_R/m_M*log(T/(J*Tr));

    return s;
}

//-----------------------------------------------------------------------------
//! specific internal energy
double FEIdealGas::SpecificInternalEnergy(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = m_Tr + tf.m_T;

    double u = SpecificFreeEnergy(mp) + T*SpecificEntropy(mp);
    
    return u;
}

//-----------------------------------------------------------------------------
//! specific gage enthalpy
double FEIdealGas::SpecificGageEnthalpy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();

    double h = SpecificInternalEnergy(mp) + m_R/m_M*(tf.m_T - fp.m_Jf*m_Tr);
    
    return h;
}

//-----------------------------------------------------------------------------
//! specific free enthalpy
double FEIdealGas::SpecificFreeEnthalpy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();

    double g = SpecificFreeEnergy(mp) + m_R/m_M*(tf.m_T - fp.m_Jf*m_Tr);
    
    return g;
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FEIdealGas::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = fp.m_Jf;
    double T = tf.m_T;
    double Tr = m_Tr;
    
    // strain-dependent contribution
    double a = m_R/m_M*(J*Tr-T+T*log(T/(J*Tr)));

    return a;
}

//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FEIdealGas::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = tf.m_T;
    
    double cp = m_C[0] + m_C[1]*T + m_C[2]*pow(T, 2) + m_C[3]*pow(T, 3) + m_C[4]*pow(T, 4);
    cp *= m_R/m_M;
    
    return cp;
}
                
//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FEIdealGas::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    double cv = IsobaricSpecificHeatCapacity(mp) - m_R/m_M;

    return cv;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FEIdealGas::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FEIdealGas::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = tf.m_T;
    
    double dcv = m_C[1] + 2*m_C[2]*T + 3*m_C[3]*pow(T, 2) + 4*m_C[4]*pow(T, 3);
    dcv *= m_R/m_M;
    
    return dcv;
}
