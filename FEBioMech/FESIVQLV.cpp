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
#include <limits>
#include <cmath>

#define sgn(x) x==0 ? 0 : x/abs(x)

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
    m_alphap = m_alpha = 1;
    m_Up = m_U = mat3dd(1);
    m_Usp = m_Us = m_U;
    m_Udp = m_Ud = m_U;
    m_R = mat3dd(1);

    // don't forget to initialize the base class
    FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
//! Update material point data.
void FESIVQLVMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    double dt = timeInfo.timeIncrement;
    m_sedp = m_sed;
    m_alphap = m_alpha;
    m_Up = m_U;
    m_Usp = m_Us;
    m_Udp = m_Ud;

    // don't forget to call the base class
    FEMaterialPointData::Update(timeInfo);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FESIVQLVMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
    ar & m_sed & m_sedp;
    ar & m_alpha & m_alphap;
    ar & m_U & m_Up;
    ar & m_Us & m_Usp;
    ar & m_Ud & m_Udp;
    ar & m_R;
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
    
    mat3dd I(1);
    mat3ds U = ep.RightStretch();
    if ((U-I).norm() <= eps)
        return s;
    mat3ds C = ep.RightCauchyGreen();
    mat3ds H = ep.RightHencky();
    mat3ds E = ep.Strain();
    double lam[3];
    vec3d u[3];
    U.eigen2(lam, u);
    mat3ds Ue[3];
    for (int i=0; i<3; ++i) Ue[i] = dyad(u[i]);
    
    // store safe copy of deformation gradient
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;
    
    // calculate alpha at current time
    double errrel = 1e-6;
    double errabs = 1e-9;
    int maxit = 100;
    int iter = 0;
    double x = pt.m_alpha; // solve for x
    bool cnvgd = false;
    bool error = false;
    mat3ds Us, Ud;
    do {
        // evaluate Us and Ud
        Us = Ud = mat3ds(0);
        mat3ds Udi(0);
        for (int i=0; i<3; ++i) {
            Us += Ue[i]*pow(lam[i],x);
            Ud += Ue[i]*pow(lam[i],1-x);
            Udi += Ue[i]/pow(lam[i],1-x);
        }
        mat3ds Cs = (Us*Us).sym();
        mat3ds Cd = (Ud*Ud).sym();
        mat3ds Cdi = (Udi*Udi).sym();
        mat3ds Es = (Cs - I)/2;
        // use Gram-Schmidt orthogonalization when evaluation U dot
        mat3ds Udot = (U*(U.dotdot(pt.m_Up)/U.dotdot(U)) - pt.m_Up)/dt;
        mat3ds Usdot = (Us*(Us.dotdot(pt.m_Usp)/Us.dotdot(Us)) - pt.m_Usp)/dt;
        mat3ds Uddot = (Ud*(Ud.dotdot(pt.m_Udp)/Ud.dotdot(Ud)) - pt.m_Udp)/dt;
        mat3ds Edot = (Udot*U).sym();
        mat3ds Esdot = (Usdot*Us).sym();
        mat3ds Eddot = (Uddot*Ud).sym();
        double alphadot = (x - pt.m_alphap)/dt;
        mat3ds Esdotc = Esdot - (Cs*H).sym()*alphadot;
        mat3ds Eddotc = Eddot + (Cd*H).sym()*alphadot;
        double Jd = Ud.det();
        mat3ds Smhat = m_Mxwl->PK2Stress(mp, Es)/(2*eta*Jd);
        double a = H.dotdot(H);
        double b = Smhat.dotdot((Cs*H).sym()) - 2*Eddotc.dotdot((H*Cdi).sym());
        double c = Smhat.dotdot((Udi*Edot*Udi).sym() - Esdotc) - (Cdi*Eddotc).dotdot(Eddotc*Cdi);
        double delta = b*b + 4*a*c;
        delta = (delta < 0) ? 0 : sqrt(delta);
        double adot = (fabs(a) <= eps) ? 0 : -(b+sgn(b)*delta)/(2*a);
        double dx = adot*dt + pt.m_alphap - x;
        x += dx;
        if (fabs(dx) <= fabs(x)*errrel) cnvgd = true;
        if (fabs(dx) <= errabs) cnvgd = true;
        if (++iter == maxit) error = true;
    } while (!cnvgd && !error);
    if (error)
        feLogWarning("SIV dashpot stretch calculation did not converge!");
    pt.m_alpha = x;
    pt.m_U = U;
    pt.m_Us = Us;
    pt.m_Ud = Ud;
    pt.m_R = Fsafe*U.inverse();

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
