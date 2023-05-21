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


#include "FELinearElasticFluid.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include "FEFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
FELinearElasticFluid::FELinearElasticFluid(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_k = 0;
    m_rhor = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FELinearElasticFluid::Init()
{
    return true;
}

//-----------------------------------------------------------------------------
// serialization
void FELinearElasticFluid::Serialize(DumpStream& ar)
{
    if (ar.IsShallow()) return;
    ar & m_k & m_rhor;
}

//-----------------------------------------------------------------------------
//! gage pressure
double FELinearElasticFluid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double p = -m_k*fp.m_ef;
    
    return p;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FELinearElasticFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return -m_k;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FELinearElasticFluid::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FELinearElasticFluid::Tangent_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FELinearElasticFluid::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FELinearElasticFluid::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FELinearElasticFluid::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double a = m_k/(2*m_rhor)*pow(fp.m_ef,2);
    return a;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FELinearElasticFluid::SpecificEntropy(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FELinearElasticFluid::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    return SpecificFreeEnergy(mp);
}

//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FELinearElasticFluid::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FELinearElasticFluid::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FELinearElasticFluid::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FELinearElasticFluid::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FELinearElasticFluid::Dilatation(const double T, const double p, double& e)
{
    e = -p/m_k;
    return true;
}
