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
#include "FESSVQLV.h"
#include "FEElasticFiberMaterial.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/DumpStream.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/matrix.h>
#include <limits>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility> // For std::pair
#include <functional> // For std::greater

//#define sgn(x) x==0 ? 0 : x/abs(x)

//-----------------------------------------------------------------------------
FESSVQLVMaterialPoint::FESSVQLVMaterialPoint(FEMaterialPointData* mp) : FEMaterialPointData(mp)
{
    m_sed = 0.0;
    m_sedp = 0.0;
    m_C = m_Cp = mat3dd(1);
    m_Cs = m_Csp = mat3dd(1);
    m_Cd = m_Cdp = mat3dd(1);
    m_G = m_Gp = mat3dd(1);
    m_sed = m_sedp = 0;
    m_Sm = m_Smp = mat3ds(0);
    m_Chatm.zero(); m_Chatmp.zero();
}

//-----------------------------------------------------------------------------
void FESSVQLVMaterialPoint::Init()
{
    m_sed = 0.0;
    m_sedp = 0.0;
    m_C = m_Cp = mat3dd(1);
    m_Cs = m_Csp = mat3dd(1);
    m_Cd = m_Cdp = mat3dd(1);
    m_G = m_Gp = mat3dd(1);
    m_sed = m_sedp = 0;
    m_Sm = m_Smp = mat3ds(0);
    m_Chatm.zero(); m_Chatmp.zero();
}

//-----------------------------------------------------------------------------
mat3ds FESSVQLVMaterialPoint::CtoH(mat3ds C)
{
    double lam2[3];
    vec3d u[3];
    C.eigen2(lam2,u);
    mat3ds H(0);
    for (int i=0; i<3; ++i) {
        H += dyad(u[i])*log(lam2[i]);
    }
    return H;
}

//-----------------------------------------------------------------------------
mat3ds FESSVQLVMaterialPoint::CtoU(mat3ds C)
{
    double lam2[3];
    vec3d u[3];
    C.eigen2(lam2,u);
    mat3ds U(0);
    for (int i=0; i<3; ++i) {
        U += dyad(u[i])*sqrt(lam2[i]);
    }
    return U;
}

//-----------------------------------------------------------------------------
mat3ds FESSVQLVMaterialPoint::HtoC(mat3ds H)
{
    double lnlam2[3];
    vec3d u[3];
    H.eigen2(lnlam2,u);
    mat3ds C(0);
    for (int i=0; i<3; ++i) {
        C += dyad(u[i])*exp(lnlam2[i]);
    }
    return C;
}

//-----------------------------------------------------------------------------
//! evaluate fourth-order tensor derivatives
tens4dmm FESSVQLVMaterialPoint::dHdC(mat3ds C)
{
    double eps = 1e-6;
    double lam2[3];
    vec3d u[3];
    C.eigen2(lam2,u);
    mat3ds Ui[3];
    for (int i=0; i<3; ++i) Ui[i] = dyad(u[i]);
    
    // check for repeated eigenvalues and reorder as needed
    int rep = 0;
    if ((fabs(lam2[1] - lam2[2]) < eps) && (fabs(lam2[0] - lam2[1]) > eps)) {
        // swap eigenvalues 2 and 0
        rep = 1;
        std::swap(lam2[0], lam2[2]);
        std::swap(u[0],u[2]);
        std::swap(Ui[0], Ui[2]);
    }
    else if ((fabs(lam2[0] - lam2[2]) < eps) && (fabs(lam2[1] - lam2[2]) > eps)) {
        // swap eigenvalues 1 and 2
        rep = 1;
        std::swap(lam2[1], lam2[2]);
        std::swap(u[1],u[2]);
        std::swap(Ui[1], Ui[2]);
    }
    else if ((fabs(lam2[0] - lam2[1]) < eps) && (fabs(lam2[1] - lam2[2]) > eps)) {
        rep = 1;
    }
    else if ((fabs(lam2[1] - lam2[2]) < eps) && (fabs(lam2[0] - lam2[1]) < eps)) {
        rep = 2;
        u[0] = vec3d(1,0,0); u[1] = vec3d(0,1,0); u[2] = vec3d(0,0,1);
        for (int i=0; i<3; ++i) Ui[i] = dyad(u[i]);
    }
    
    tens4dmm result, UixUi;
    result.zero(); UixUi.zero();
    for (int i=0; i<3; ++i) {
        UixUi += dyad1mm(Ui[i]);
        result += dyad1mm(Ui[i]/lam2[i]);
    }
    
    tens4dmm dUidC[3];
    if (rep == 0) {
        dUidC[0] = dyad4mm(Ui[0], Ui[1])/(lam2[0]-lam2[1]) + dyad4mm(Ui[0], Ui[2])/(lam2[0]-lam2[2]);
        dUidC[1] = dyad4mm(Ui[1], Ui[2])/(lam2[1]-lam2[2]) + dyad4mm(Ui[1], Ui[0])/(lam2[1]-lam2[0]);
        dUidC[2] = dyad4mm(Ui[2], Ui[0])/(lam2[2]-lam2[0]) + dyad4mm(Ui[2], Ui[1])/(lam2[2]-lam2[1]);
        
        for (int i=0; i<3; ++i) result += dUidC[i]*log(lam2[i]);
    }
    else if (rep == 1) {
        dUidC[2] = dyad4mm(Ui[2], Ui[0])/(lam2[2]-lam2[0]) + dyad4mm(Ui[2], Ui[1])/(lam2[2]-lam2[1]);
        result += dUidC[2]*(log(lam2[2]) - log(lam2[0]));
    }
    
    return result;
}

tens4dmm FESSVQLVMaterialPoint::dCdH(mat3ds C)
{
    double eps = 1e-6;
    double lam2[3];
    vec3d u[3];
    C.eigen2(lam2,u);
    mat3ds Ui[3];
    for (int i=0; i<3; ++i) Ui[i] = dyad(u[i]);
    
    // check for repeated eigenvalues and reorder as needed
    int rep = 0;
    if ((fabs(lam2[1] - lam2[2]) < eps) && (fabs(lam2[0] - lam2[1]) > eps)) {
        // swap eigenvalues 2 and 0
        rep = 1;
        std::swap(lam2[0], lam2[2]);
        std::swap(Ui[0], Ui[2]);
    }
    else if ((fabs(lam2[0] - lam2[2]) < eps) && (fabs(lam2[1] - lam2[2]) > eps)) {
        // swap eigenvalues 1 and 2
        rep = 1;
        std::swap(lam2[1], lam2[2]);
        std::swap(Ui[1], Ui[2]);
    }
    else if ((fabs(lam2[0] - lam2[1]) < eps) && (fabs(lam2[1] - lam2[2]) > eps)) {
        rep = 1;
    }
    else if ((fabs(lam2[1] - lam2[2]) < eps) && (fabs(lam2[0] - lam2[1]) < eps)) {
        rep = 2;
        u[0] = vec3d(1,0,0); u[1] = vec3d(0,1,0); u[2] = vec3d(0,0,1);
        for (int i=0; i<3; ++i) Ui[i] = dyad(u[i]);
    }
    
    double lnlam2[3];
    for (int i=0; i<3; ++i) lnlam2[i] = log(lam2[i]);
    tens4dmm result, UixUi, IoI;
    IoI = dyad4mm(mat3dd(1));
    result.zero(); UixUi.zero();
    for (int i=0; i<3; ++i) {
        tens4dmm tmp = dyad1mm(Ui[i]);
        UixUi += tmp;
        result += tmp*lam2[i];
    }
    
    tens4dmm dUidH[3];
    if (rep == 0) {
        dUidH[0] = dyad4mm(Ui[0], Ui[1])/(lnlam2[0]-lnlam2[1]) + dyad4mm(Ui[0], Ui[2])/(lnlam2[0]-lnlam2[2]);
        dUidH[1] = dyad4mm(Ui[1], Ui[2])/(lnlam2[1]-lnlam2[2]) + dyad4mm(Ui[1], Ui[0])/(lnlam2[1]-lnlam2[0]);
        dUidH[2] = dyad4mm(Ui[2], Ui[0])/(lnlam2[2]-lnlam2[0]) + dyad4mm(Ui[2], Ui[1])/(lnlam2[2]-lnlam2[1]);
        
        for (int i=0; i<3; ++i) result += dUidH[i]*lam2[i];
    }
    else if (rep == 1) {
        dUidH[2] = dyad4mm(Ui[2], Ui[0])/(lnlam2[2]-lnlam2[0]) + dyad4mm(Ui[2], Ui[1])/(lnlam2[2]-lnlam2[1]);
        result += dUidH[2]*(lam2[2]-lam2[0]);
    }
    
    return result;
}

//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPointData* FESSVQLVMaterialPoint::Copy()
{
    FESSVQLVMaterialPoint* pt = new FESSVQLVMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
//! p p material point data.
void FESSVQLVMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    double dt = timeInfo.timeIncrement;
    m_sedp = m_sed;
    m_Cp = m_C;
    m_Csp = m_Cs;
    m_Cdp = m_Cd;
    m_Gp = m_G;
    m_Smp = m_Sm;
    m_Chatmp = m_Chatm;

    // don't forget to call the base class
    FEMaterialPointData::Update(timeInfo);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FESSVQLVMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
    ar & m_sed & m_sedp;
    ar & m_C & m_Cp;
    ar & m_Cs & m_Csp;
    ar & m_Cd & m_Cdp;
    ar & m_G & m_Gp;
    ar & m_Sm & m_Smp;
    ar & m_Chatm & m_Chatmp;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESSVQLV, FEElasticMaterial)

    // material parameters
    ADD_PARAMETER(m_eta  , "eta"  )->setLongName("Dashpot shear viscosity")->setUnits(UNIT_VISCOSITY);
    ADD_PARAMETER(m_kappa, "kappa")->setLongName("Dashpot bulk viscosity" )->setUnits(UNIT_VISCOSITY);

    // define the material properties
    ADD_PROPERTY(m_Base, "parallel");
    ADD_PROPERTY(m_Mxwl, "Maxwell");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESSVQLV::FESSVQLV(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_Base = nullptr;
    m_Mxwl = nullptr;
    m_eta = 0.0;
    m_kappa = 0.0;
    m_secant_tangent = true;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPointData* FESSVQLV::CreateMaterialPointData()
{
    return new FESSVQLVMaterialPoint(m_Base->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! update specialize material point data
void FESSVQLV::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FESSVQLV::Stress(FEMaterialPoint& mp)
{
    FETimeInfo& timePoint = GetFEModel()->GetTime();
    double t = timePoint.currentTime;
    double dt = timePoint.timeIncrement;
    double tp = t - dt;
    
    double errrel = 1e-6;
    double errabs = 1e-9;
    int itmax = 100;
    
    double eta = m_eta(mp);
    double kappa = m_kappa(mp);

    // get the elastic material point data
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    // get the viscoelastic point data
    FESSVQLVMaterialPoint& pt = *mp.ExtractData<FESSVQLVMaterialPoint>();

    mat3dd I(1);
    tens4dmm IoI = dyad4mm(I);
    
    // evaluate base stress
    pt.m_C = ep.RightCauchyGreen();
    mat3ds H = pt.CtoH(pt.m_C);
    if (H.norm() < errrel) return mat3ds(0);
    mat3ds E = (pt.m_C-I)/2;
    mat3ds S = m_Base->PK2Stress(mp, E);

    bool convgd = false;
    bool maxed = false;
    int it = 0;
    mat3ds Hsp = pt.CtoH(pt.m_Csp);
    tens4dmm dCsdHs;
    // initialize Cs to previoust time point
    pt.m_Cs = pt.m_C;
    double Jd;
    do {
        // Evaluate Cd
        mat3ds Hs = pt.CtoH(pt.m_Cs);
        mat3ds Hd = H - Hs;
        pt.m_Cd = pt.HtoC(Hd);
        Jd = sqrt(pt.m_Cd.det());
        dCsdHs = pt.dCdH(pt.m_Cs);
        tens4dmm dHsdCs = pt.dHdC(pt.m_Cs);
        tens4dmm tmpdot = ddot(dCsdHs,dHsdCs);
        tens4d dCsdHs4d(dCsdHs);
        tens4d dCsdHsi = dCsdHs4d.inverse();
        mat3ds Gd = ddot(dCsdHs,pt.dHdC(pt.m_Cd)).dot(pt.m_Cd);
        pt.m_G = ddot(dCsdHs,pt.dHdC(pt.m_C)).dot(pt.m_C);
        double tmp = (kappa > 0) ? 1./(3*kappa) : 0;
        pt.m_Chatm = (dyad4mm(Gd)/2/eta+dyad1mm(Gd)*(tmp-1./2/eta)/3)*(2./ep.m_J);
        mat3ds Es = (pt.m_Cs-I)/2;
        pt.m_Sm = m_Mxwl->PK2Stress(mp, Es);
        mat3ds dCs = -pt.m_Cs;
        pt.m_Cs = pt.m_Csp + pt.m_G - pt.m_Gp-(pt.m_Chatm.dot(pt.m_Sm) + pt.m_Chatmp.dot(pt.m_Smp))*dt/2;
        dCs += pt.m_Cs;
        if (fabs(dCs.dotdot(pt.m_Cs)) <= errrel) convgd = true;
        if (dCs.norm() <= errabs) convgd = true;
        if (++it > itmax) { convgd = true; maxed = true; }
    } while (!convgd);
    
    if (maxed) feLogWarning("SSV-QLV iterations did not converge!\n");
    
    S += ddot(pt.dHdC(pt.m_C),dCsdHs).dot(pt.m_Sm);
    mat3ds s = (ep.m_F*(S*(ep.m_F).transpose())/ep.m_J).sym();

    // evaluate residual dissipation
//    pt.m_rd = ?;
    
    // return the total Cauchy stress,
    return s;
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FESSVQLV::Tangent(FEMaterialPoint& mp)
{
    // we should not come here since we are evaluating the tangent numerically
    assert(false);
    tens4ds c; c.zero();
    return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density function
double FESSVQLV::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    if (dt == 0) return 0;
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the viscoelastic point data
    FESSVQLVMaterialPoint& pt = *mp.ExtractData<FESSVQLVMaterialPoint>();

    // Calculate the new elastic Cauchy stress
    double sed = m_Base->StrainEnergyDensity(mp);

    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;

    // evaluate Fs and Js to calculate stress in Maxwell spring
    mat3ds Us = pt.CtoU(pt.m_Cs);
    ep.m_F = Us; // Fs
    ep.m_J = Us.det();  // Js
    sed += m_Mxwl->StrainEnergyDensity(mp);

    ep.m_F = Fsafe;
    ep.m_J = Jsafe;

    // return the total Cauchy stress,
    return sed;
}
