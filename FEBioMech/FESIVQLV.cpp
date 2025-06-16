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
    m_Udp = m_Ud = mat3dd(1);
    m_R = mat3dd(1);
    m_Udotp = m_Udot = mat3ds(0);
    m_Omegap = m_Omega = mat3da(vec3d(0,0,0));
    m_SEw = 0;
    m_alphap = m_alpha = 1;

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
    m_Udp = m_Ud;
    m_Udotp = m_Udot;
    m_Omegap = m_Omega;
    m_alphap = m_alpha;

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
    ar & m_Ud & m_Udp;
    ar & m_R;
    ar & m_Udot & m_Udotp;
    ar & m_Omega & m_Omegap;
    ar & m_SEw;
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
    double errrel = 1e-6;
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
    
    double lam[3], H[3];
    vec3d u[3];
    mat3ds e(0), edot(0);
    mat3dd I(1);
    mat3ds U = ep.RightStretch();   // U at current time
    mat3ds C = ep.RightCauchyGreen();   // C at current time
    mat3ds E = ep.Strain();         // E at current time
    mat3ds Udot = (U - pt.m_Up)/dt; // Udot over interval
    if ((U-I).norm() <= errabs) {
        pt.m_Usp.eigen2(lam,u);
        for (int i=0; i<3; ++i) {
            u[i] = (I + pt.m_Omegap)*u[i];
            u[i].Normalize();
            lam[i] = 1;
            H[i] = 0;
        }
    }
    else {
        U.eigen2(lam,u);
        int imax = 0;
        for (int i=0; i<3; ++i) {
            H[i] = log(lam[i]);
            if (fabs(H[i]) > fabs(H[imax])) imax = i;
        }
        // reorder eigenvalues and eigenvectors
        double ltmp = lam[2], htmp = H[2];
        vec3d utmp = u[2];
        lam[2] = lam[imax]; H[2] = H[imax];
        u[2] = u[imax];
        lam[imax] = ltmp; H[imax] = htmp;
        u[imax] = utmp;
    }

    mat3ds Ue[3];
    double lamdot[3], lams[3], lamsp[3];
    double Hdot[3];
    double Hs[3], Hsp[3];
    mat3ds Udot0(0);
    for (int i=0; i<3; ++i) {
        u[i].Normalize();
        Ue[i] = dyad(u[i]);
        lamdot[i] = Udot.dotdot(Ue[i]);
        Hdot[i] = lamdot[i]/lam[i];
        Udot0 += Ue[i]*lamdot[i];
        lams[i] = pt.m_Us.dotdot(Ue[i]);
        lamsp[i] = pt.m_Usp.dotdot(Ue[i]);
        Hs[i] = log(lams[i]);
        Hsp[i] = log(lamsp[i]);
    }
    mat3ds OUUO = Udot - Udot0;
    vec3d omtmp(0,0,0);
    if (OUUO.norm() > errabs) {
        double dlam = lam[0]-lam[1];
        if (fabs(dlam) > errabs) omtmp.z = u[1]*(OUUO*u[0])/dlam;
        dlam = lam[1]-lam[2];
        if (fabs(dlam) > errabs) omtmp.x = u[2]*(OUUO*u[1])/dlam;
        dlam = lam[2]-lam[0];
        if (fabs(dlam) > errabs) omtmp.y = u[0]*(OUUO*u[2])/dlam;
    }
    vec3d omega = u[0]*omtmp.x + u[1]*omtmp.y + u[2]*omtmp.z;
    mat3da Omega(omega);
    
    // store safe copy of deformation gradient
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;
    
    mat3ds Us = pt.m_Us;
    mat3ds Ud;
    int iter = 0;
    bool cnvgd = false;
    bool error = false;
    double SEw = 0;
    do {
        vector<double> f(3,0.);
        matrix df(3,3);
        double lamd[3];
        Us = mat3ds(0);
        Ud = mat3ds(0);
        mat3ds Udi(0);
        for (int i=0; i<3; ++i) {
            lams[i] = exp(Hs[i]);
            lamd[i] = lam[i]/lams[i];
            Us += Ue[i]*lams[i];
            Ud += Ue[i]*lamd[i];
            Udi += Ue[i]/lamd[i];
        }
        double Jd = lamd[0]*lamd[1]*lamd[2];
        mat3ds Es = ((Us*Us).sym()-I)/2;
        mat3ds Smhat = m_Mxwl->PK2Stress(mp, Es)/(2*eta*Jd);
        mat3ds Ewdot = 2*((Udi*Omega*E*Udi).sym() - (Omega*Es).sym());
        SEw = Smhat.dotdot(Ewdot);
        double a[3], b[3], c[3];
        double dx = 0;
        if (fabs(H[2]) > errabs) {
            for (int i=0; i<3; ++i) {
                a[i] = dt*pow(H[2],2)*Hdot[i]+H[i]*H[2]*Hsp[2];
                b[i] = dt*(H[2]*Hdot[i]-H[i]*Hdot[2])+H[i]*H[2];
                c[i] = a[i] - dt*pow(H[2]*lams[i],2)*Smhat.dotdot(Ue[i]);
            }
            double a3 = 0;
            double b3 = 0;
            double c3 = -pow(dt,2)*pow(H[2],4)*SEw;
            for (int i=0; i<3; ++i) {
                a3 += pow(b[i],2);
                b3 += (c[i]+a[i])*b[i];
                c3 += c[i]*a[i];
            }
            double delta = b3*b3 - 4*a3*c3;
            if (delta < 0)
                delta = 0;
            else
                delta = sqrt(delta);
            double root1 = (b3 + delta)/(2*a3);
            double root2 = (b3 - delta)/(2*a3);
            double root;
            if (b3 < 0)
                root = root1;
            else
                root = root2;
            dx = root - Hs[2];
            Hs[2] += dx;
            Hs[0] = H[0]/H[2]*Hs[2];
            Hs[1] = H[1]/H[2]*Hs[2];
        }
        else {
            dx = Hs[2] - (Hsp[2] + dt*Hdot[2]);
            Hs[0] = Hsp[0] + dt*Hdot[0];
            Hs[1] = Hsp[1] + dt*Hdot[1];
            Hs[2] = Hsp[2] + dt*Hdot[2];
        }
        double x = Hs[2];
        if (dx <= fabs(x)*errrel) cnvgd = true;
        if (++iter == maxit) error = true;
        if (error)
            bool pause = true;
    } while (!cnvgd && !error);
    if (error)
        feLogWarning("SIV dashpot stretch calculation did not converge!");
    for (int i=0; i<3; ++i) lams[i] = exp(Hs[i]);
    pt.m_Us = Ue[0]*lams[0] + Ue[1]*lams[1] + Ue[2]*lams[2];
    pt.m_U = U;
    pt.m_Omega = Omega;
    pt.m_Udot = Udot;
    pt.m_SEw = SEw;
    pt.m_R = Fsafe*U.inverse();
    pt.m_Ud = (U*pt.m_Us.inverse()).sym();
    pt.m_alpha = 0;
    if (fabs(H[2]) > errabs) pt.m_alpha = Hs[2]/H[2];

    // evaluate Fs and Js to calculate stress in Maxwell spring
    ep.m_F = pt.m_R*pt.m_Us;     // Fs
    ep.m_J = lams[0]*lams[1]*lams[2]; // Js
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
