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
#include "FESIVQLV.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/DumpStream.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/matrix.h>
#include <limits>
#include <cmath>

//#define sgn(x) x==0 ? 0 : x/abs(x)

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESIVQLV, FEElasticMaterial)

    // material parameters
    ADD_PARAMETER(m_eta, "eta")->setUnits(UNIT_VISCOSITY);

    // define the material properties
    ADD_PROPERTY(m_Base, "parallel");
    ADD_PROPERTY(m_Mxwl, "Maxwell");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESIVQLVMaterialPoint::FESIVQLVMaterialPoint(FEMaterialPointData* mp) : FEMaterialPointData(mp)
{
    m_sed = 0.0;
    m_sedp = 0.0;
}

//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPointData* FESIVQLVMaterialPoint::Copy()
{
    FESIVQLVMaterialPoint* pt = new FESIVQLVMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FESIVQLVMaterialPoint::Init()
{
    // intialize data to zero
    m_sedp = m_sed = 0.0;
    m_Up = m_U = mat3dd(1);
    m_Usp = m_Us = mat3dd(1);
    m_Ed = m_Edp = mat3ds(0);
    m_R = mat3dd(1);
    m_Udotp = m_Udot = mat3ds(0);

    // don't forget to initialize the base class
    FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
//! Update material point data.
void FESIVQLVMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    double dt = timeInfo.timeIncrement;
    m_sedp = m_sed;
    m_Up = m_U;
    m_Usp = m_Us;
    m_Edp = m_Ed;
    m_Udotp = m_Udot;

    // don't forget to call the base class
    FEMaterialPointData::Update(timeInfo);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FESIVQLVMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
    ar & m_sed & m_sedp;
    ar & m_U & m_Up;
    ar & m_Us & m_Usp;
    ar & m_Ed & m_Edp;
    ar & m_R;
    ar & m_Udot & m_Udotp;
}

//-----------------------------------------------------------------------------
//! constructor
FESIVQLV::FESIVQLV(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_eta = 0;
    m_Base = nullptr;
    m_Mxwl = nullptr;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPointData* FESIVQLV::CreateMaterialPointData()
{
    return new FESIVQLVMaterialPoint(m_Base->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! update specialize material point data
void FESIVQLV::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FESIVQLV::Stress(FEMaterialPoint& mp)
{
    const double eps = 10*std::numeric_limits<double>::epsilon();
    double errrel = 1e-3;
    double errabs = 1e-9;
    int maxit = 100;

    FETimeInfo& tp = GetFEModel()->GetTime();
    double dt = tp.timeIncrement;
    if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
    // Calculate the base Cauchy stress
    mat3ds s = m_Base->Stress(mp);
    
    double eta = m_eta(mp);
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the viscoelastic point data
    FESIVQLVMaterialPoint& pt = *mp.ExtractData<FESIVQLVMaterialPoint>();
    
    double lam[3];
    vec3d u[3];
    mat3ds e(0), edot(0);
    mat3dd I(1);
    mat3ds U = ep.RightStretch();   // U at current time
    mat3ds Udot = (U - pt.m_Up)/dt; // Udot over interval
    U.eigen2(lam, u);               // eigenvalues & eigenvectors at current time

    mat3ds Ue[3];
    double lamdot[3];
    mat3da Omega(vec3d(0,0,0));
    mat3ds Udot0(0);
    for (int i=0; i<3; ++i) {
        Ue[i] = dyad(u[i]);
        lamdot[i] = Udot.dotdot(Ue[i]);
        Udot0 += Ue[i]*lamdot[i];
    }

    mat3ds OUUO = Udot - Udot0;
    double omg[3];
    double dlam;
    dlam = lam[0] - lam[1]; omg[2] = (dlam != 0) ? u[1]*(OUUO*u[0])/dlam : 0;
    dlam = lam[1] - lam[2]; omg[0] = (dlam != 0) ? u[2]*(OUUO*u[1])/dlam : 0;
    dlam = lam[2] - lam[0]; omg[1] = (dlam != 0) ? u[0]*(OUUO*u[2])/dlam : 0;
    vec3d omega (omg[0], omg[1], omg[2]);
    Omega = mat3da(omega);
    
    // store safe copy of deformation gradient
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;
    
    // calculate Ed0-dot
    // and use it to calculate Ed
    mat3ds Ed = pt.m_Ed;
    int iter = 0;
    bool cnvgd = false;
    bool error = false;
    mat3ds Us;
    double lamd[3], lams[3];
    do {
        Us = mat3ds(0);
        // evaluate the eigenvalues of Ud and Us
        for (int i=0; i<3; ++i) {
            lamd[i] = sqrt(1+2*Ed.dotdot(Ue[i]));
            lams[i] = lam[i]/lamd[i];
            Us += Ue[i]*lams[i];
        }
        double Jdm = lamd[0]*lamd[1]*lamd[2];
        mat3ds Es = ((Us*Us).sym()-I)/2;
        mat3ds Smhat = m_Mxwl->PK2Stress(mp, Es)/(2*eta*Jdm);
        mat3ds Ed0dot(0), Edp(0);
        for (int i=0; i<3; ++i) {
            Ed0dot += Ue[i]*(Smhat.dotdot(Ue[i])*lam[i]);
            Edp += Ue[i]*pt.m_Edp.dotdot(Ue[i]);
        }
        mat3ds dEd = Edp + Ed0dot*dt - Ed;
        Ed += dEd;
        double dEdn = dEd.norm();
        if (dEdn <= Ed.norm()*errrel) cnvgd = true;
        if (dEdn <= errabs) cnvgd = true;
        if (++iter == maxit) error = true;
    } while (!cnvgd && !error);
    if (error)
        feLogWarning("SIV dashpot stretch calculation did not converge!");
    pt.m_Ed = Ed;
    pt.m_U = U;
    pt.m_Udot = Udot;
    pt.m_R = Fsafe*U.inverse();
    pt.m_Us = Us;
    pt.m_alpha = 0;
    int imax = 0;
    double lmax = fabs(log(lam[imax]));
    for (int i=0; i<3; ++i) {
        double l = fabs(log(lam[i]));
        if (l > lmax) { imax = i; lmax = l;}
    }
    if (lmax > 0) pt.m_alpha = log(lams[imax])/log(lam[imax]);

    // evaluate Fs and Js to calculate stress in Maxwell spring
    ep.m_F = pt.m_R*pt.m_Us;     // Fs
    ep.m_J = pt.m_Us.det(); // Js
    s += m_Mxwl->Stress(mp)*(ep.m_J/Jsafe);
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;

    // return the total Cauchy stress,
    return s;
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FESIVQLV::Tangent(FEMaterialPoint& mp)
{
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
}

//-----------------------------------------------------------------------------
//! Strain energy density function
double FESIVQLV::StrainEnergyDensity(FEMaterialPoint& mp)
{
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
}
