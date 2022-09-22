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
#include "FEThermoFluidMaterialPoint.h"
#include "FEThermoFluid.h"

//-----------------------------------------------------------------------------
//! specific internal energy
double FEElasticFluid::SpecificInternalEnergy(FEMaterialPoint& mp)
{
    double Tr = GetGlobalConstant("T");

    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double T = tf.m_T + Tr;

    double u = SpecificFreeEnergy(mp) + T*SpecificEntropy(mp);
    
    return u;
}

//-----------------------------------------------------------------------------
//! specific gage enthalpy
double FEElasticFluid::SpecificGageEnthalpy(FEMaterialPoint& mp)
{
    FEThermoFluid* pMat = dynamic_cast<FEThermoFluid*>(GetParent());
    
    double h = SpecificInternalEnergy(mp) + Pressure(mp)/pMat->Density(mp);
    
    return h;
}

//-----------------------------------------------------------------------------
//! specific free enthalpy
double FEElasticFluid::SpecificFreeEnthalpy(FEMaterialPoint& mp)
{
    FEThermoFluid* pMat = dynamic_cast<FEThermoFluid*>(GetParent());
    
    double g = SpecificFreeEnergy(mp) + Pressure(mp)/pMat->Density(mp);
    
    return g;
}

//-----------------------------------------------------------------------------
//! pressure from state variables
double FEElasticFluid::Pressure(const double ef, const double T)
{
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = ef;
    ft->m_T = T;
    FEMaterialPoint tmp(ft);
    double p = Pressure(tmp);
    delete ft;
    return p;
}

//-----------------------------------------------------------------------------
//! dp/dJ
double FEElasticFluid::Tangent_Strain(const double ef, const double T)
{
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = ef;
    ft->m_T = T;
    FEMaterialPoint tmp(ft);
    double dpJ = Tangent_Strain(tmp);
    delete ft;
    return dpJ;
}

//-----------------------------------------------------------------------------
//! dp/dT
double FEElasticFluid::Tangent_Temperature(const double ef, const double T)
{
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = ef;
    ft->m_T = T;
    FEMaterialPoint tmp(ft);
    double dpT = Tangent_Temperature(tmp);
    delete ft;
    return dpT;
}

//-----------------------------------------------------------------------------
//! d2p/dJ2
double FEElasticFluid::Tangent_Strain_Strain(const double ef, const double T)
{
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = ef;
    ft->m_T = T;
    FEMaterialPoint tmp(ft);
    double dpJ2 = Tangent_Strain_Strain(tmp);
    delete ft;
    return dpJ2;
}

//-----------------------------------------------------------------------------
//! d2p/dJdT
double FEElasticFluid::Tangent_Strain_Temperature(const double ef, const double T)
{
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = ef;
    ft->m_T = T;
    FEMaterialPoint tmp(ft);
    double dpJT = Tangent_Strain_Temperature(tmp);
    delete ft;
    return dpJT;
}

//-----------------------------------------------------------------------------
//! d2p/dT2
double FEElasticFluid::Tangent_Temperature_Temperature(const double ef, const double T)
{
    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = ef;
    ft->m_T = T;
    FEMaterialPoint tmp(ft);
    double dpT2 = Tangent_Temperature_Temperature(tmp);
    delete ft;
    return dpT2;
}
