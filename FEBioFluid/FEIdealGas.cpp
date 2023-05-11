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
    ADD_PARAMETER(m_M   , FE_RANGE_GREATER(0.0), "M")->setUnits("M/n")->setLongName("molar mass");
    ADD_PARAMETER(m_ar  , "ar")->setLongName("normalized referential specific free energy");    // ar normalized by R.Tr/M
    ADD_PARAMETER(m_sr  , "sr")->setLongName("normalized referential specific entropy");        // sr normalized by R/M
    ADD_PROPERTY (m_ao  , "ao")->SetLongName("normalized specific free energy circle");         // a-circle normalized by R.Tr/M
    ADD_PROPERTY (m_cp  , "cp")->SetLongName("normalized isobaric specific heat capacity");     // cp normalized by R/M

END_FECORE_CLASS();

FEIdealGas::FEIdealGas(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_R = m_Pr = m_Tr = m_ar = m_sr = 0;
    m_ao = nullptr;
    m_cp = nullptr;
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
void FEIdealGas::Serialize(DumpStream& ar)
{
    FEElasticFluid::Serialize(ar);
    if (ar.IsShallow()) return;
    
    ar & m_R & m_Pr & m_Tr;
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
    double That = T/m_Tr;
    double scl = m_R*m_Tr/m_M;
    
    // referential free energy
    double a = m_ar - m_sr*(That-1);
    
    // add a_circle
    a += m_ao->value(That);
    
    // add strain-dependent contribution
    a += J+That*(log(That/J)-1);

    return a*scl;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FEIdealGas::SpecificEntropy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double scl = m_R/m_M;

    // referential entropy
    double s = m_sr;
    
    // add s_circle
    s -= m_ao->derive(That);
    
    // add strain-dependent contribution
    s += -log(That/J);

    return s*scl;
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FEIdealGas::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double T = tf.m_T;
    double That = T/m_Tr;
    double scl = m_R*m_Tr/m_M;

    // strain-dependent contribution
    double a = J+That*(log(That/J)-1);

    return a*scl;
}

//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FEIdealGas::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = tf.m_T + m_Tr;
    double That = T/m_Tr;
    double scl = m_R/m_M;

    double cp = m_cp->value(That);
    
    return cp*scl;
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
    double That = T/m_Tr;
    double scl = m_R/(m_M*m_Tr);

    double dcv = m_cp->derive(That);
    
    return dcv*scl;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FEIdealGas::Dilatation(const double T, const double p, double& e)
{
    double J = (T+m_Tr)/m_Tr/(1+p/m_Pr);
    e = J - 1;
    return true;
}
