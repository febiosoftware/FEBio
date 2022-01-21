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


#include "FEIdealGas.h"
#include <FECore/log.h>
#include "FEFluidMaterialPoint.h"
#include "FEThermoFluid.h"
#include "FEThermoFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEIdealGas, FEElasticFluid)

    // material parameters
    ADD_PARAMETER(m_M   , FE_RANGE_GREATER(0.0), "M");
    ADD_PARAMETER(m_ar  , "ar");
    ADD_PARAMETER(m_sr  , "sr");
    ADD_PROPERTY (m_ao  , "ao");
    ADD_PROPERTY (m_cp  , "cp");

END_FECORE_CLASS();

FEIdealGas::FEIdealGas(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_R = m_Pr = m_Tr = m_ar = m_sr = 0;
    m_ao = m_cp = nullptr;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEIdealGas::Init()
{
    m_R  = GetGlobalConstant("R");
    m_Tr = GetGlobalConstant("T");
    m_Pr = GetGlobalConstant("P");
    
    if (m_R  <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section");    return false; }
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    if (m_Pr == 0) {
        FEThermoFluid* pMat = dynamic_cast<FEThermoFluid*>(GetParent());
        double rhor = pMat->ReferentialDensity();
        m_Pr = m_R*m_Tr/m_M*rhor;
        feLogWarning("The referential absolute pressure P is calculated internally as %g\n",m_Pr);
    }
    
    m_ao->Init();
    m_cp->Init();
    
    return true;
}

//-----------------------------------------------------------------------------
//! gage pressure
double FEIdealGas::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = m_Tr + tf.m_T;
    double Jf = 1 + fp.m_ef;

    double p = m_Pr*(T/(Jf*m_Tr) - 1);

    return p;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FEIdealGas::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = m_Tr + tf.m_T;
    double Jf = 1 + fp.m_ef;

    double dp = -m_Pr*T/(Jf*Jf*m_Tr);

    return dp;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FEIdealGas::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = m_Tr + tf.m_T;
    double Jf = 1 + fp.m_ef;

    double d2p = 2*m_Pr*T/(pow(Jf,3)*m_Tr);

    return d2p;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FEIdealGas::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double Jf = 1 + fp.m_ef;

    double dp = m_Pr/(Jf*m_Tr);

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
    double Jf = 1 + fp.m_ef;

    double d2p = -m_Pr/(Jf*Jf*m_Tr);

    return d2p;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FEIdealGas::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double Tr = m_Tr;
    
    // referential free energy
    double a = m_ar - m_sr*(T-Tr);
    
    // add a_circle
    a += m_ao->value(T);
    
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
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double Tr = m_Tr;
    
    // referential entropy
    double s = m_sr;
    
    // add s_circle
    s -= m_ao->derive(T);
    
    // add strain-dependent contribution
    s += -m_R/m_M*log(T/(J*Tr));

    return s;
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FEIdealGas::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
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
    double T = tf.m_T + m_Tr;
    
    double cp = m_cp->value(T);
    
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
    double T = tf.m_T + m_Tr;
    
    double dcv = m_cp->derive(T);
    
    return dcv;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FEIdealGas::Dilatation(const double T, const double p, const double c, double& e)
{
    double J = (T+m_Tr)/m_Tr/(1+p/m_Pr);
    e = J - 1;
    return true;
}
