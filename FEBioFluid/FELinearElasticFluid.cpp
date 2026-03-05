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
#include "FEFluid.h"

//-----------------------------------------------------------------------------
FELinearElasticFluid::FELinearElasticFluid(FEModel* pfem) : FEElasticFluid(pfem)
{
}

//-----------------------------------------------------------------------------
//! initialization
bool FELinearElasticFluid::Init()
{
    return true;
}

//-----------------------------------------------------------------------------
//! gauge pressure
double FELinearElasticFluid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
	double k = m_pFluid->m_k;
    double p = -k*fp.m_ef;
    
    return p;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FELinearElasticFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return -m_pFluid->m_k;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FELinearElasticFluid::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FELinearElasticFluid::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
	double k = m_pFluid->m_k;
	double rhor = m_pFluid->m_rhor;

    double a = k/(2*rhor)*pow(fp.m_ef,2);
    return a;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FELinearElasticFluid::Dilatation(const double T, const double p, double& e)
{
    e = -p/m_pFluid->m_k;
    return true;
}
