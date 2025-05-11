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
#include "FESIVQLVUC.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/DumpStream.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <limits>
#include <cmath>

#define sgn(x) x==0 ? 0 : x/abs(x)

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESIVQLVUC, FEUncoupledMaterial)

    // material parameters
    ADD_PARAMETER(m_eta, "eta")->setUnits(UNIT_VISCOSITY);

    // define the material properties
    ADD_PROPERTY(m_Base, "parallel");
    ADD_PROPERTY(m_Mxwl, "Maxwell");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESIVQLVUC::FESIVQLVUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
    m_eta = 0;
    m_Base = nullptr;
    m_Mxwl = nullptr;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPointData* FESIVQLVUC::CreateMaterialPointData()
{
    return new FESIVQLVMaterialPoint(m_Base->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! update specialize material point data
void FESIVQLVUC::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FESIVQLVUC::DevStress(FEMaterialPoint& mp)
{
    const double eps = 10*std::numeric_limits<double>::epsilon();
    
    FETimeInfo& tp = GetFEModel()->GetTime();
    double dt = tp.timeIncrement;
    if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
    
    // Calculate the base Cauchy stress
    mat3ds s = m_Base->DevStress(mp);
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the viscoelastic point data
    FESIVQLVMaterialPoint& pt = *mp.ExtractData<FESIVQLVMaterialPoint>();
    
    // store safe copy of deformation gradient
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;
    
    // calculate principal dashpot stretch
    // terms are accumulated in s
    double errrel = 1e-6;
    double errabs = 1e-15;
    int maxit = 100;
    int iter = 0;
    double x = pt.m_lamdp; // solve for x
    bool cnvgd = false;
    bool error = false;
    mat3ds U = ep.RightStretch();
    mat3ds C = ep.RightCauchyGreen();
    mat3dd I(1);
    if (tp.currentTime > 0.5)
        bool pause = true;
    do {
        mat3dd Ud(x);
        mat3ds Us = U/x;
        mat3ds Es = (C/(x*x)-I)/2;
        double Jd = Ud.det();
        mat3ds Smhat = m_Mxwl->PK2Stress(mp, Es)/(6*m_eta);
        double c = Smhat.dotdot(C);
        double g = x - pt.m_lamdp - dt*pow(x,-4)*c;
        double dg = 1 + dt*pow(x,-5)*4*c;
        double dx = -g/dg;
        x += dx;
        if (fabs(dx) <= fabs(x)*errrel) cnvgd = true;
        if (fabs(g) <= errabs) cnvgd = true;
        if (++iter == maxit) { error = true; }
    } while (!cnvgd && !error);
    if (error)
        feLogWarning("SIV dashpot stretch calculation did not converge!");
    pt.m_lamd = x;
    pt.m_U = U;
    pt.m_Us = U/x;
    pt.m_Ud = mat3dd(x);
    pt.m_R = Fsafe*U.inverse();
    
    // evaluate Fs and Js to calculate stress in Maxwell spring
    ep.m_F = pt.m_R*pt.m_Us;     // Fs
    ep.m_J = pt.m_Us.det(); // Js
    s += m_Mxwl->Stress(mp)*(ep.m_J/Jsafe);
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;
    
    // return the total Cauchy stress,
    return s.dev();
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FESIVQLVUC::DevTangent(FEMaterialPoint& mp)
{
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // store safe copy of deformation gradient
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;
    
    // get the viscoelastic point data
    FESIVQLVMaterialPoint& pt = *mp.ExtractData<FESIVQLVMaterialPoint>();
    
    // Calculate the new elastic Cauchy stress
    tens4ds c = m_Base->DevTangent(mp);
    mat3ds s = m_Base->DevStress(mp);

    // evaluate Fs and Js to calculate stress in Maxwell spring
    ep.m_F = pt.m_R*pt.m_Us; // Fs
    ep.m_J = pt.m_Us.det();  // Js
    c += m_Mxwl->Tangent(mp)*(ep.m_J/Jsafe);
    s += m_Mxwl->Stress(mp)*(ep.m_J/Jsafe);

    ep.m_F = Fsafe;
    ep.m_J = Jsafe;
    
    // This is the final value of the elasticity tensor
    mat3dd I(1);
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    
    c += - 1./3.*(ddots(c,IxI) - IxI*(c.tr()/3.))
    + 2./3.*((I4-IxI/3.)*s.tr()-dyad1s(s.dev(),I));
    
    // return the total elastic tangent,
    return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density function
double FESIVQLVUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    if (dt == 0) return 0;
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the viscoelastic point data
    FESIVQLVMaterialPoint& pt = *mp.ExtractData<FESIVQLVMaterialPoint>();
    
    // Calculate the new elastic Cauchy stress
    double sed = m_Base->DevStrainEnergyDensity(mp);
    
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;
    
    // evaluate Fs and Js to calculate stress in Maxwell spring
    ep.m_F = pt.m_R*pt.m_Us; // Fs
    ep.m_J = pt.m_Us.det();  // Js
    sed += m_Mxwl->StrainEnergyDensity(mp);
    
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;
    
    // return the total Cauchy stress,
    return sed;
}
