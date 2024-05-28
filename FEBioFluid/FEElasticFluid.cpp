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
//! tangent of pressure with respect to strain J
double FEElasticFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    double d = 1e-6;
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();

    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = pf.m_ef+d;
    ft->m_T = tf.m_T;
    FEMaterialPoint tmp(ft);
    double pp = Pressure(tmp);
    fp->m_ef = pf.m_ef-d;
    double pm = Pressure(tmp);
    delete ft;
    double dpJ = (pp - pm)/(2*d);
    return dpJ;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FEElasticFluid::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    double d = 1e-6;
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();

    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = pf.m_ef+d;
    ft->m_T = tf.m_T;
    FEMaterialPoint tmp(ft);
    double pp = Pressure(tmp);
    fp->m_ef = pf.m_ef-d;
    double pm = Pressure(tmp);
    fp->m_ef = pf.m_ef;
    double p = Pressure(tmp);
    delete ft;
    double dpJ2 = (pp - 2*p + pm)/(d*d);
    return dpJ2;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FEElasticFluid::Tangent_Temperature(FEMaterialPoint& mp)
{
    double Tr = GetGlobalConstant("T");
    double d = 1e-6*Tr;
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();

    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = pf.m_ef;
    ft->m_T = tf.m_T+d;
    FEMaterialPoint tmp(ft);
    double pp = Pressure(tmp);
    ft->m_T = tf.m_T-d;
    double pm = Pressure(tmp);
    delete ft;
    double dpT = (pp - pm)/(2*d);
    return dpT;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FEElasticFluid::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    double Tr = GetGlobalConstant("T");
    double d = 1e-6*Tr;
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();

    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = pf.m_ef;
    ft->m_T = tf.m_T+d;
    FEMaterialPoint tmp(ft);
    double pp = Pressure(tmp);
    ft->m_T = tf.m_T-d;
    double pm = Pressure(tmp);
    ft->m_T = tf.m_T;
    double p = Pressure(tmp);
    delete ft;
    double dpT2 = (pp - 2*p + pm)/(d*d);
    return dpT2;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FEElasticFluid::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    double Tr = GetGlobalConstant("T");
    double dJ = 1e-6;
    double dT = 1e-6*Tr;
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();

    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = pf.m_ef+dJ;
    ft->m_T = tf.m_T+dT;
    FEMaterialPoint tmp(ft);
    double ppp = Pressure(tmp);
    fp->m_ef = pf.m_ef-dJ;
    double ppm = Pressure(tmp);
    ft->m_T = tf.m_T-dT;
    double pmm = Pressure(tmp);
    fp->m_ef = pf.m_ef+dJ;
    double pmp = Pressure(tmp);

    delete ft;
    double dpTJ = (ppp - ppm - pmp + pmm)/(4*dT*dJ);
    return dpTJ;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FEElasticFluid::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    double d = 1e-6;
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();

    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = pf.m_ef+d;
    ft->m_T = tf.m_T;
    FEMaterialPoint tmp(ft);
    double cvp = IsochoricSpecificHeatCapacity(tmp);
    fp->m_ef = pf.m_ef-d;
    double cvm = IsochoricSpecificHeatCapacity(tmp);
    delete ft;
    double dcvJ = (cvp - cvm)/(2*d);
    return dcvJ;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FEElasticFluid::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    double Tr = GetGlobalConstant("T");
    double d = 1e-6*Tr;
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();

    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = pf.m_ef;
    ft->m_T = tf.m_T+d;
    FEMaterialPoint tmp(ft);
    double cvp = IsochoricSpecificHeatCapacity(tmp);
    ft->m_T = tf.m_T-d;
    double cvm = IsochoricSpecificHeatCapacity(tmp);
    delete ft;
    double dcvT = (cvp - cvm)/(2*d);
    return dcvT;
}

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
//! specific gauge enthalpy
double FEElasticFluid::SpecificGaugeEnthalpy(FEMaterialPoint& mp)
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
