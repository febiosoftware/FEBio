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


#include "FEIdealLiquid.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include "FEFluidMaterialPoint.h"
#include "FEThermoFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEIdealLiquid, FEElasticFluid)

    // material parameters
    ADD_PARAMETER(m_k   , FE_RANGE_GREATER(0.0)         , "k" );
    ADD_PARAMETER(m_cv  , FE_RANGE_GREATER_OR_EQUAL(0.0), "cv");
    ADD_PARAMETER(m_beta, "beta");
    ADD_PARAMETER(m_ar  , "ar");
    ADD_PARAMETER(m_sr  , "sr");

END_FECORE_CLASS();

FEIdealLiquid::FEIdealLiquid(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_R = m_Pr = m_Tr = 0;
    m_k = 0;
    m_beta = 0;
    m_cv = 0;
    m_ar = m_sr = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEIdealLiquid::Init()
{
    m_R  = GetFEModel()->GetGlobalConstant("R");
    m_Tr = GetFEModel()->GetGlobalConstant("T");
    m_Pr = GetFEModel()->GetGlobalConstant("P");
    
    if (m_R  <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section");    return false; }
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    if (m_Pr <= 0) { feLogError("A positive referential absolute pressure P must be defined in Globals section"); return false; }
    
    m_pMat = dynamic_cast<FEThermoFluid*>(GetParent());
    m_rhor = m_pMat->m_rhor;
    
    return true;
}

//-----------------------------------------------------------------------------
//! gage pressure
double FEIdealLiquid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double p = m_k*(1-fp.m_Jf) + m_beta*(tf.m_T);

    return p;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FEIdealLiquid::Tangent_Strain(FEMaterialPoint& mp)
{
    return -m_k;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FEIdealLiquid::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FEIdealLiquid::Tangent_Temperature(FEMaterialPoint& mp)
{
    return m_beta;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FEIdealLiquid::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FEIdealLiquid::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FEIdealLiquid::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = fp.m_Jf;
    double T = tf.m_T + m_Tr;
    double Tr = m_Tr;
    
    double a = m_k/(2*m_rhor)*pow(J - 1, 2) - m_beta/m_rhor*(T - Tr)*(J - 1)
    + (m_cv - m_sr)*(T - Tr) - m_cv*T*log(T/Tr)+ m_ar;
    
    return a;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FEIdealLiquid::SpecificEntropy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = fp.m_Jf;
    double T = tf.m_T + m_Tr;
    double Tr = m_Tr;
    
    // referential entropy
    double s = m_beta/m_rhor*(J - 1) + m_cv*log(T/Tr) + m_sr;
    
    return s;
}

//-----------------------------------------------------------------------------
//! specific internal energy
double FEIdealLiquid::SpecificInternalEnergy(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = tf.m_T + m_Tr;

    double u = SpecificFreeEnergy(mp) + T*SpecificEntropy(mp);
    
    return u;
}

//-----------------------------------------------------------------------------
//! specific gage enthalpy
double FEIdealLiquid::SpecificGageEnthalpy(FEMaterialPoint& mp)
{
    double h = SpecificInternalEnergy(mp) + Pressure(mp)/m_pMat->Density(mp);
    
    return h;
}

//-----------------------------------------------------------------------------
//! specific free enthalpy
double FEIdealLiquid::SpecificFreeEnthalpy(FEMaterialPoint& mp)
{
    double g = SpecificFreeEnergy(mp) + Pressure(mp)/m_pMat->Density(mp);
    
    return g;
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FEIdealLiquid::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = fp.m_Jf;
    
    // strain-dependent contribution
    double a = m_k/(2*m_rhor)*pow(J - 1, 2) - m_beta/m_rhor*tf.m_T*(J - 1);

    return a;
}

//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FEIdealLiquid::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    return m_cv;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FEIdealLiquid::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FEIdealLiquid::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    return 0;
}
