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



#include "FEIdealGasIsentropic.h"
#include "FEFluidMaterialPoint.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEIdealGasIsentropic, FEElasticFluid)
	ADD_PARAMETER(m_gamma, FE_RANGE_GREATER(0.0), "gamma")->setLongName("specific heat capacity ratio");
	ADD_PARAMETER(m_M    , FE_RANGE_GREATER(0.0), "M"    )->setLongName("molar mass")->setUnits(UNIT_MOLAR_MASS);
    ADD_PARAMETER(m_cv0  , FE_RANGE_GREATER(0.0), "cv0"  )->setLongName("isochoric specific heat capacity")->setUnits(UNIT_SPECIFIC_ENTROPY);
END_FECORE_CLASS();

//============================================================================
// FEIdealGasIsentropic
//============================================================================

//-----------------------------------------------------------------------------
//! FEIdealGasIsentropic constructor

FEIdealGasIsentropic::FEIdealGasIsentropic(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_rhor = 0;
    m_k = 0;
    m_gamma = 0;
    m_M = 0;
    m_cv0 = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEIdealGasIsentropic::Init() 
{
    m_R  = GetGlobalConstant("R");
    m_Tr = GetGlobalConstant("T");
    m_Pr = GetGlobalConstant("P");
    
	if (m_R  <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section"); return false; }
	if (m_Tr <= 0) { feLogError("A positive ambient absolute temperature T must be defined in Globals section"); return false; }
	if (m_Pr <= 0) { feLogError("A positive ambient absolute pressure P must be defined in Globals section"); return false; }

    if (m_rhor == 0) m_rhor = m_M*m_Pr/(m_R*m_Tr);
    if (m_k == 0) m_k = m_Pr;
    
    return true;
}

//-----------------------------------------------------------------------------
void FEIdealGasIsentropic::Serialize(DumpStream& ar)
{
    FEElasticFluid::Serialize(ar);
    if (ar.IsShallow()) return;
    
    ar & m_R & m_Pr & m_Tr & m_rhor;
}

//-----------------------------------------------------------------------------
//! elastic pressure
double FEIdealGasIsentropic::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double p =m_k*(pow(1+fp.m_ef,-m_gamma)-1);
    return p;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FEIdealGasIsentropic::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double dpJ = -m_k*m_gamma*pow(1+fp.m_ef,-1-m_gamma);
    return dpJ;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FEIdealGasIsentropic::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double d2pJ = m_k*m_gamma*(1+m_gamma)*pow(1+fp.m_ef,-2-m_gamma);
    return d2pJ;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FEIdealGasIsentropic::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double dpT = m_k/(1+fp.m_ef)/m_Tr;
    return dpT;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FEIdealGasIsentropic::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FEIdealGasIsentropic::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double dpTJ = -m_k/pow(1+fp.m_ef,2)/m_Tr;
    return 0;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FEIdealGasIsentropic::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double a = m_Tr*(m_R/m_M*fp.m_ef + m_cv0*(pow(1+fp.m_ef,1-m_gamma)-1));
    return a;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FEIdealGasIsentropic::SpecificEntropy(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FEIdealGasIsentropic::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    return SpecificFreeEnergy(mp);
}

//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FEIdealGasIsentropic::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    return m_cv0 + m_R/m_M;
}

//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FEIdealGasIsentropic::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    return m_cv0;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FEIdealGasIsentropic::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FEIdealGasIsentropic::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FEIdealGasIsentropic::Dilatation(const double T, const double p, double& e)
{
    e = pow(1+p/m_k,-1/m_gamma)-1.0;
    return true;
}
