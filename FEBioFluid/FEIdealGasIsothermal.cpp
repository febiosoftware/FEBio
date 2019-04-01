/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
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
    m_R  = GetFEModel()->GetGlobalConstant("R");
    m_Tr = GetFEModel()->GetGlobalConstant("T");
    m_pr = GetFEModel()->GetGlobalConstant("p");
    
    if (m_R  <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section");    return false; }
    if (m_pr <= 0) { feLogError("A positive ambient absolute pressure p must be defined in Globals section"); return false; }
    
    m_rhor = m_M*m_pr/(m_R*m_Tr);
    
    return true;
}

//-----------------------------------------------------------------------------
//! elastic pressure from dilatation
double FEIdealGasIsothermal::Pressure(const double e)
{
    double J = 1 + e;
    return m_pr*(1./J - 1);
}

//-----------------------------------------------------------------------------
//! tangent of elastic pressure with respect to strain J
double FEIdealGasIsothermal::Tangent_Pressure_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = fp.m_Jf;
    double dp = -m_pr/J;
    return dp;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of elastic pressure with respect to strain J
double FEIdealGasIsothermal::Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = fp.m_Jf;
    double d2p = 2*m_pr/(J*J);
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
    double J = fp.m_Jf;
    double sed = m_pr*(J-1-log(J));
    return sed;
}

//-----------------------------------------------------------------------------
//! invert pressure-dilatation relation
double FEIdealGasIsothermal::Dilatation(const double p)
{
    double J = m_pr/(p+m_pr);
    return J - 1;
}

