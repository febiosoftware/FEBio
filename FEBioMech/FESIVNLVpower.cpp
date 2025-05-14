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
#include "FESIVNLVpower.h"
#include "FEUncoupledMaterial.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/DumpStream.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <limits>
#include <cmath>

#define sgn(x) x==0 ? 0 : x/abs(x)

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESIVNLVpower, FEElasticMaterial)

    // material parameters
    ADD_PARAMETER(m_e0, "eta0")->setUnits(UNIT_VISCOSITY);
    ADD_PARAMETER(m_e1, "eta1")->setUnits(UNIT_VISCOSITY);
    ADD_PARAMETER(m_E0, "E0")->setUnits(UNIT_NONE);
    ADD_PARAMETER(m_a , "alpha")->setUnits(UNIT_NONE);
    ADD_PARAMETER(m_bdash, "dashpot");

    // define the material properties
    ADD_PROPERTY(m_Base, "parallel");
    ADD_PROPERTY(m_Mxwl, "Maxwell");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESIVNLVpower::FESIVNLVpower(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_e0 = m_e1 = 0;
    m_E0 = 1;
    m_a = 1;
    m_bdash = true;
    m_Base = nullptr;
    m_Mxwl = nullptr;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPointData* FESIVNLVpower::CreateMaterialPointData()
{
    return new FESIVQLVMaterialPoint(m_Base->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! update specialize material point data
void FESIVNLVpower::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FESIVNLVpower::Stress(FEMaterialPoint& mp) { return mat3ds(0); }
/*{
    const double eps = 10*std::numeric_limits<double>::epsilon();
    
    FETimeInfo& tp = GetFEModel()->GetTime();
    double dt = tp.timeIncrement;
    if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
    
    // Calculate the base Cauchy stress
    mat3ds s = m_Base->Stress(mp);
    
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
    do {
        mat3dd Ud(x);
        mat3ds Us = U/x;
        mat3ds Es = (C/(x*x)-I)/2;
        mat3ds Ed = mat3dd((x*x-1)/2);
        double Jd = Ud.det();
        double Enorm = m_bdash ? Ed.norm() : Es.norm();
        // power-law relation for viscosity zeta as a function of dashpot/spring strain Enorm
        double eta = m_e0 + m_e1*pow(Enorm/m_E0,m_a);
        mat3ds Smhat = m_Mxwl->PK2Stress(mp, Es)/(6*eta);
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
    return s;

}*/

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FESIVNLVpower::Tangent(FEMaterialPoint& mp) { return tens4ds(0.); }
/*{
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // store safe copy of deformation gradient
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;
    
    // get the viscoelastic point data
    FESIVQLVMaterialPoint& pt = *mp.ExtractData<FESIVQLVMaterialPoint>();
    
    // Calculate the new elastic Cauchy stress
    tens4ds c = m_Base->Tangent(mp);
    
    // evaluate Fs and Js to calculate stress in Maxwell spring
    ep.m_F = pt.m_R*pt.m_Us; // Fs
    ep.m_J = pt.m_Us.det();  // Js
    c += m_Mxwl->Tangent(mp)*(ep.m_J/Jsafe);
    
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;
    
    // return the total elastic tangent,
    return c;
}*/

//-----------------------------------------------------------------------------
//! Strain energy density function
double FESIVNLVpower::StrainEnergyDensity(FEMaterialPoint& mp) { return 0; }
/*{
    double dt = GetFEModel()->GetTime().timeIncrement;
    if (dt == 0) return 0;
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the viscoelastic point data
    FESIVQLVMaterialPoint& pt = *mp.ExtractData<FESIVQLVMaterialPoint>();
    
    // Calculate the new elastic Cauchy stress
    double sed = m_Base->StrainEnergyDensity(mp);
    
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
}*/
