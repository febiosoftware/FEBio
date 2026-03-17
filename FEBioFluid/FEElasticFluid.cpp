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
#include "FEElasticFluid.h"
#include "FEFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
// By default, the absolute temperature is the ambient temperature for isothermal fluids
double FEElasticFluid::Temperature(FEMaterialPoint& mp)
{
    return GetGlobalConstant("T");
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FEElasticFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    double d = 1e-6;
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();

	FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    fp->m_ef = pf.m_ef+d;
    FEMaterialPoint tmp(fp);
    double pp = Pressure(tmp);
    fp->m_ef = pf.m_ef-d;
    double pm = Pressure(tmp);
    double dpJ = (pp - pm)/(2*d);
    return dpJ;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FEElasticFluid::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    double d = 1e-6;
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();

    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    fp->m_ef = pf.m_ef+d;
    FEMaterialPoint tmp(fp);
    double dpp = Tangent_Strain(tmp);
    fp->m_ef = pf.m_ef-d;
    double dpm = Tangent_Strain(tmp);
    double dpJ2 = (dpp - dpm)/(2*d);
    return dpJ2;
}

//-----------------------------------------------------------------------------
//! pressure from state variables
double FEElasticFluid::Pressure(const double ef, const double T)
{
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    fp->m_ef = ef;
    FEMaterialPoint tmp(fp);
    double p = Pressure(tmp);
    return p;
}

//-----------------------------------------------------------------------------
//! dp/dJ
double FEElasticFluid::Tangent_Strain(const double ef, const double T)
{
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    fp->m_ef = ef;
    FEMaterialPoint tmp(fp);
    double dpJ = Tangent_Strain(tmp);
    return dpJ;
}
