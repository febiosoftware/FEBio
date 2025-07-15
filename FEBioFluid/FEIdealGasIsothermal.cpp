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



#include "stdafx.h"
#include "FEIdealGasIsothermal.h"
#include "FEFluidMaterialPoint.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEIdealGasIsothermal, FEElasticFluid)
ADD_PARAMETER(m_M    , FE_RANGE_GREATER(0.0), "M"    )->setLongName("molar mass")->setUnits(UNIT_MOLAR_MASS);
ADD_PARAMETER(m_cv0  , FE_RANGE_GREATER(0.0), "cv0"  )->setLongName("isochoric specific heat capacity")->setUnits(UNIT_SPECIFIC_ENTROPY);
END_FECORE_CLASS();

//============================================================================
// FEIdealGasIsothermal
//============================================================================

//-----------------------------------------------------------------------------
//! FEIdealGasIsentropic constructor

FEIdealGasIsothermal::FEIdealGasIsothermal(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_rhor = 0;
    m_M = 0;
    m_cv0 = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEIdealGasIsothermal::Init()
{
    m_R  = GetGlobalConstant("R");
    m_Tr = GetGlobalConstant("T");
    m_Pr = GetGlobalConstant("P");
    
    if (m_R  <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section"); return false; }
    if (m_Tr <= 0) { feLogError("A positive ambient absolute temperature T must be defined in Globals section"); return false; }
    if (m_Pr <= 0) { feLogError("A positive ambient absolute pressure P must be defined in Globals section"); return false; }
    
    if (m_rhor == 0) m_rhor = m_M*m_Pr/(m_R*m_Tr);
    
    return true;
}

//-----------------------------------------------------------------------------
void FEIdealGasIsothermal::Serialize(DumpStream& ar)
{
    FEElasticFluid::Serialize(ar);
    if (ar.IsShallow()) return;
    
    ar & m_R & m_Pr & m_Tr & m_rhor;
}

//-----------------------------------------------------------------------------
//! elastic pressure
double FEIdealGasIsothermal::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;
    double phat = 1/J-1;
    return phat*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FEIdealGasIsothermal::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;
    double dphatJ = -1/pow(J,2);
    return dphatJ*m_Pr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FEIdealGasIsothermal::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;
    double d2phatJ = 2/pow(J,3);
    return d2phatJ*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FEIdealGasIsothermal::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;
    double dphatT = 1/J;
    return dphatT*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FEIdealGasIsothermal::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FEIdealGasIsothermal::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;
    double d2phatJT = -1/pow(J,2);
    return d2phatJT*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FEIdealGasIsothermal::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;
    double ahat = J + (log(1/J)-1);
    return ahat*m_Pr/m_rhor;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FEIdealGasIsothermal::SpecificEntropy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;
    double shat = -log(1/J);
    return shat*m_R/m_M;
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FEIdealGasIsothermal::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;
    double what = J + (log(1/J)-1);
    return what*m_Pr/m_rhor;
}

//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FEIdealGasIsothermal::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    return m_cv0 + m_R/m_M;
}

//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FEIdealGasIsothermal::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    return m_cv0;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FEIdealGasIsothermal::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FEIdealGasIsothermal::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FEIdealGasIsothermal::Dilatation(const double T, const double p, double& e)
{
    e = 1/(1+p/m_Pr)-1;
    return true;
}
