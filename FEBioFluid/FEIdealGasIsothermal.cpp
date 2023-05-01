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
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEIdealGasIsothermal, FEFluid)
	ADD_PARAMETER(m_M    , FE_RANGE_GREATER(0.0), "M"    );
END_FECORE_CLASS();

//============================================================================
// FEIdealGasIsothermal
//============================================================================

//-----------------------------------------------------------------------------
//! FEIdealGasIsothermal constructor

FEIdealGasIsothermal::FEIdealGasIsothermal(FEModel* pfem) : FEFluid(pfem)
{
    m_rhor = 0;
    m_k = 0;
    m_M = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEIdealGasIsothermal::Init()
{
    m_R  = GetGlobalConstant("R");
    m_Tr = GetGlobalConstant("T");
    m_Pr = GetGlobalConstant("P");
    
    if (m_R  <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section");    return false; }
    if (m_Pr <= 0) { feLogError("A positive ambient absolute pressure P must be defined in Globals section"); return false; }
    
    m_rhor = m_M*m_Pr/(m_R*m_Tr);
    
    return true;
}

//-----------------------------------------------------------------------------
void FEIdealGasIsothermal::Serialize(DumpStream& ar)
{
    FEFluid::Serialize(ar);
    if (ar.IsShallow()) return;
    
    ar & m_R & m_Pr & m_Tr & m_rhor;
}

//-----------------------------------------------------------------------------
//! elastic pressure
double FEIdealGasIsothermal::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    return Pressure(fp.m_ef,0);
}

//-----------------------------------------------------------------------------
//! elastic pressure from dilatation
double FEIdealGasIsothermal::Pressure(const double e, const double T)
{
    double J = 1 + e;
    return m_Pr*(1./J - 1);
}

//-----------------------------------------------------------------------------
//! tangent of elastic pressure with respect to strain J
double FEIdealGasIsothermal::Tangent_Pressure_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;
    double dp = -m_Pr/J;
    return dp;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of elastic pressure with respect to strain J
double FEIdealGasIsothermal::Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;
    double d2p = 2*m_Pr/(J*J);
    return d2p;
}

//-----------------------------------------------------------------------------
//! evaluate temperature
double FEIdealGasIsothermal::Temperature(FEMaterialPoint& mp)
{
    return m_Tr;
}

//-----------------------------------------------------------------------------
//! calculate free energy density (per reference volume)
double FEIdealGasIsothermal::StrainEnergyDensity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;
    double sed = m_Pr*(J-1-log(J));
    return sed;
}

//-----------------------------------------------------------------------------
//! invert effective pressure-dilatation relation
bool FEIdealGasIsothermal::Dilatation(const double T, const double p, double& e)
{
    double J = m_Pr/(p+m_Pr);
    e = J - 1;
    return true;
}

