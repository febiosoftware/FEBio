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
    m_E = m_Ep = mat3dd(0);
    m_Es = m_Esp = mat3dd(0);
    m_G = m_Gp = mat3dd(0);
    m_sed = m_sedp = 0;
}

//-----------------------------------------------------------------------------
void FESSVQLVMaterialPoint::Init()
{
    m_E = m_Ep = mat3dd(0);
    m_Es = m_Esp = mat3dd(0);
    m_G = m_Gp = mat3dd(0);
    m_sed = m_sedp = 0;
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
    m_Ep = m_E;
    m_Esp = m_Es;
    m_Gp = m_G;

    // don't forget to call the base class
    FEMaterialPointData::Update(timeInfo);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FESSVQLVMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
    ar & m_sed & m_sedp;
    ar & m_E & m_Ep;
    ar & m_Es & m_Esp;
    ar & m_G & m_Gp;
}

//-----------------------------------------------------------------------------
//! Evaluate right stretch tensor from right Cauchy-Green tensor
mat3ds FESSVQLVMaterialPoint::CtoU(mat3ds C) {
    double lam[3];
    vec3d u[3];
    C.eigen2(lam,u);
    mat3ds U; U.zero();
    for (int i=0; i<3; ++i) {
        lam[i] = sqrt(lam[i]);
        U += dyads(u[i])*lam[i];
    }
    return U;
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
    
    double errrel = 1e-9;
    double errabs = 1e-12;
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
    pt.m_E = ep.Strain();
    if (pt.m_E.norm() < errrel) return mat3ds(0);
    mat3ds S = m_Base->PK2Stress(mp, pt.m_E);

    bool convgd = false;
    bool maxed = false;
    int it = 0;
    mat3ds Sm; Sm.zero();
    // initialize Es to previoust time point
    pt.m_Es = (pt.m_Esp.norm() > errrel) ? pt.m_Esp : pt.m_E;
    do {
        // Evaluate Cd
        mat3ds Ed = pt.m_E - pt.m_Es;
        mat3ds Cd = I + Ed*2;
        // Evaluate G
        double tmpk = (kappa > 0) ? 1./(3*kappa) : 0;
        double tmpe = (eta > 0) ? 1./(2*eta) : 0;
        tens4dmm Hdi = (dyad4s(Cd)*tmpe + dyad1s(Cd)*((tmpk - tmpe)/3))/ep.m_J;
        // Evaluate G
        Sm = m_Mxwl->PK2Stress(mp,pt.m_Es);
        pt.m_G = Hdi.dot(Sm);
        mat3ds dEs = -pt.m_Es;
        pt.m_Es = pt.m_Esp + pt.m_E - pt.m_Ep - (pt.m_G + pt.m_Gp)*(dt/2);
        dEs += pt.m_Es;
        if (fabs(dEs.dotdot(pt.m_Es)) <= errrel) convgd = true;
        if (dEs.norm() <= errabs) convgd = true;
        if (++it > itmax) { convgd = true; maxed = true; }
        // make sure we get at least two iterations
        if (it < 2) { convgd = false; maxed = false; }
    } while (!convgd);
    
    if (maxed) feLogWarning("SSV-QLV iterations did not converge!\n");
    
    S += Sm;
    mat3ds s = (ep.m_F*(S*(ep.m_F).transpose())/ep.m_J).sym();

    // evaluate residual dissipation
    mat3ds Eddot = (pt.m_E - pt.m_Ep - pt.m_Es + pt.m_Esp)/dt;
    pt.m_rd = Sm.dotdot(Eddot)/ep.m_J;
    
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
    mat3ds Cs = mat3dd(1) + pt.m_Es*2;
    ep.m_F = pt.CtoU(Cs); // Fs
    ep.m_J = ep.m_F.det();  // Js
    sed += m_Mxwl->StrainEnergyDensity(mp);

    ep.m_F = Fsafe;
    ep.m_J = Jsafe;

    // return the total Cauchy stress,
    return sed;
}
