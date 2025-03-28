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
#include "FESIVViscoelastic.h"
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
BEGIN_FECORE_CLASS(FESIVViscoelastic, FEElasticMaterial)

    // material parameters
    ADD_PARAMETER(m_t[0], "t1")->setUnits(UNIT_TIME);
    ADD_PARAMETER(m_t[1], "t2")->setUnits(UNIT_TIME);
    ADD_PARAMETER(m_t[2], "t3")->setUnits(UNIT_TIME);
    ADD_PARAMETER(m_t[3], "t4")->setUnits(UNIT_TIME);
    ADD_PARAMETER(m_t[4], "t5")->setUnits(UNIT_TIME);
    ADD_PARAMETER(m_t[5], "t6")->setUnits(UNIT_TIME);
    ADD_PARAMETER(m_g0  , "g0");
    ADD_PARAMETER(m_g[0], "g1");
    ADD_PARAMETER(m_g[1], "g2");
    ADD_PARAMETER(m_g[2], "g3");
    ADD_PARAMETER(m_g[3], "g4");
    ADD_PARAMETER(m_g[4], "g5");
    ADD_PARAMETER(m_g[5], "g6");
    ADD_PARAMETER(m_ttype, FE_RANGE_CLOSED(0,2), "trigger");

    // define the material properties
    ADD_PROPERTY(m_Base, "parallel");
    ADD_PROPERTY(m_Mxwl, "Maxwell");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESIVViscoelasticMaterialPoint::FESIVViscoelasticMaterialPoint(FEMaterialPointData* mp) : FEMaterialPointData(mp)
{
    m_sed = 0.0;
    m_sedp = 0.0;
}

//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPointData* FESIVViscoelasticMaterialPoint::Copy()
{
    FESIVViscoelasticMaterialPoint* pt = new FESIVViscoelasticMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FESIVViscoelasticMaterialPoint::Init()
{
    // intialize data to zero
    m_sed = 0.0;
    m_sedp = 0.0;
    m_R = mat3dd(1);
    m_up[0] = m_u[0] = vec3d(1,0,0);
    m_up[1] = m_u[1] = vec3d(0,1,0);
    m_up[2] = m_u[2] = vec3d(0,0,1);
    for (int i=0; i<3; ++i) {
        m_U[i] = dyad(m_u[i]);
        m_lam[i] = 1.0;
    }
    for (int i=0; i<MAX_TERMS; ++i) {
        m_lam3d[i] = m_lam3dp[i] = 1.0;
        m_HAmr[i] = 0;
    }

    // don't forget to initialize the base class
    FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
//! Update material point data.
void FESIVViscoelasticMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    double dt = timeInfo.timeIncrement;
    m_sedp = m_sed;
    
    // copy previous data
    for (int i=0; i<MAX_TERMS; ++i) {
        m_lam3dp[i] = m_lam3d[i];
    }
    for (int i=0; i<3; ++i) {
        m_up[i] = m_u[i];
    }

    // don't forget to call the base class
    FEMaterialPointData::Update(timeInfo);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FESIVViscoelasticMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
    ar & m_sed & m_sedp;
    ar & m_lam3d & m_lam3dp;
    ar & m_U;
    ar & m_u & m_up;
    ar & m_lam;
    ar & m_R;
    ar & m_HAmr;
}

//-----------------------------------------------------------------------------
//! constructor
FESIVViscoelastic::FESIVViscoelastic(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_g0 = 1;
    for (int i=0; i<MAX_TERMS; ++i)
    {
        m_t[i] = 1;
        m_g[i] = 0;
    }
    m_ttype = 0;

    m_Base = nullptr;
    m_Mxwl = nullptr;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPointData* FESIVViscoelastic::CreateMaterialPointData()
{
    return new FESIVViscoelasticMaterialPoint(m_Base->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! update specialize material point data
void FESIVViscoelastic::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    // get the viscosities in the reference configuration
    FESIVViscoelasticMaterialPoint& pt = *mp.ExtractData<FESIVViscoelasticMaterialPoint>();
    mat3dd I(1.);
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3ds Ui = ep.RightStretchInverse();
    mat3d R = ep.m_F*Ui;
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;

    // set current deformation to pure rotation
    ep.m_F = R;
    ep.m_J = 1;
    mat3ds U3 = pt.m_U[2];
    tens4ds Cr = m_Mxwl->Tangent(mp);
    double HA3 = (Cr.dot(U3)).dotdot(U3);
    // restore current deformation
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;
    
    for (int i=0; i<MAX_TERMS; ++i) {
        if (m_g[i] > 0) {
            pt.m_HAmr[i] = HA3;
        }
    }
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FESIVViscoelastic::Stress(FEMaterialPoint& mp)
{
    const double eps = 10*std::numeric_limits<double>::epsilon();
    
    FETimeInfo& tp = GetFEModel()->GetTime();
    double dt = tp.timeIncrement;
    if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    

    // Calculate the base Cauchy stress
    mat3ds s = m_Base->Stress(mp)*m_g0;
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the viscoelastic point data
    FESIVViscoelasticMaterialPoint& pt = *mp.ExtractData<FESIVViscoelasticMaterialPoint>();
    
    // evaluate total strain
    mat3dd I(1);
    mat3ds C = ep.RightCauchyGreen();
    double lam2[3], lam[3];
    vec3d u[3];
    mat3ds U[3];
    C.eigen2(lam2, u);
    lam[0] = sqrt(lam2[0]); lam[1] = sqrt(lam2[1]); lam[2] = sqrt(lam2[2]);
    U[0] = dyad(u[0]); U[1] = dyad(u[1]); U[2] = dyad(u[2]);
    // evaluate rotation tensor from inverse of stretch
    mat3ds Ui = U[0]/lam[0] + U[1]/lam[1] + U[2]/lam[2];
    pt.m_R = ep.m_F*Ui;
    // identify direction of loading from largest magnitude strains
    double Ev[3];
    Ev[0] = (lam2[0]-1)/2; Ev[1] = (lam2[1]-1)/2; Ev[2] = (lam2[2]-1)/2;
    bool reorder = false;
    if ((fabs(Ev[0]) > eps) || (fabs(Ev[1]) > eps) || (fabs(Ev[2]) > eps)) {
        // reorder principal stretches
        // let dir[2] be the principal direction of largest strain magnitude
        if ((fabs(Ev[1]) > fabs(Ev[0])) && (fabs(Ev[1]) > fabs(Ev[2]))){
            pt.m_lam[2] = lam[1]; pt.m_lam[0] = lam[2]; pt.m_lam[1] = lam[0];
            pt.m_U[2] = U[1]; pt.m_U[0] = U[2]; pt.m_U[1] = U[0];
            pt.m_u[2] = u[1]; pt.m_u[0] = u[2]; pt.m_u[1] = u[0];
        }
        else if ((fabs(Ev[2]) > fabs(Ev[0])) && (fabs(Ev[2]) > fabs(Ev[1]))) {
            pt.m_lam[2] = lam[2]; pt.m_lam[0] = lam[0]; pt.m_lam[1] = lam[1];
            pt.m_U[2] = U[2]; pt.m_U[0] = U[0]; pt.m_U[1] = U[1];
            pt.m_u[2] = u[2]; pt.m_u[0] = u[0]; pt.m_u[1] = u[1];
        }
        else {
            pt.m_lam[2] = lam[0]; pt.m_lam[0] = lam[1]; pt.m_lam[1] = lam[2];
            pt.m_U[2] = U[0]; pt.m_U[0] = U[1]; pt.m_U[1] = U[2];
            pt.m_u[2] = u[0]; pt.m_u[0] = u[1]; pt.m_u[1] = u[2];
        }
        reorder = true;
    }
    // store safe copy of deformation gradient
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;
    
    // calculate principal dashpot stretch
    // terms are accumulated in s
    double errrel = 1e-3;
    double errabs = 1e-15;
    int maxit = 100;
    int iter = 0;
    if (reorder) {
        vec3d u3dot = (pt.m_u[2] - pt.m_up[2])/dt;
        double u3dot2 = u3dot*u3dot;
        //    if ((tp.timeStep == 0) && (tp.currentIteration == 0))
        //        u3dot2 = 0;
        mat3ds U3dot = dyads(u3dot, pt.m_u[2]);
        for (int i=0; i<MAX_TERMS; ++i)
        {
            // only solve for dashpot stretch if there is a non-zero time constant and gamma associated with this term
            if ((m_t[i] > 0) && (m_g[i] > 0)) {
                double lams[3]; // Maxwell spring principal stretches
                bool cnvgd = false;
                bool error = false;
                double f = 0;
                double fp = 0;
                do {
                    fp = f;
                    double x = pt.m_lam3d[i];
                    lams[0] = pt.m_lam[0]; lams[1] = pt.m_lam[1]; lams[2] = pt.m_lam[2]/x;
                    mat3ds Us = pt.m_U[0]*lams[0] + pt.m_U[1]*lams[1] + pt.m_U[2]*lams[2];
                    mat3ds Es = ((Us*Us).sym() - I)/2;
                    ep.m_F = pt.m_R*Us;
                    ep.m_J = Us.det();
                    mat3ds Smhat = m_Mxwl->PK2Stress(mp,Es)/pt.m_HAmr[i];
                    double c1 = Smhat.dotdot(pt.m_U[2])/m_t[i];
                    double c2 = Smhat.dotdot(U3dot)/2/m_t[i];
                    double c3 = Smhat.dotdot((U3dot*C + C*U3dot).sym())/2/m_t[i];
                    double c4 = Smhat.dotdot((pt.m_U[2]*U3dot*C + C*U3dot*pt.m_U[2]).sym())/2/m_t[i];
                    double lam3 = pt.m_lam[2];
                    double lam3dp = pt.m_lam3dp[i];
                    double delta = pow(c1,2)*pow(lam3,4) -
                    4*c2*pow(lam3,2)*pow(-1 + x,2)*pow(x,6)*(1 + x) -
                    (-1 + x)*pow(x,7)*(4*c4*(-1 + x) - 4*c3*x + u3dot2*(-1 + x)*x*pow(1 + x,2));
                    if (delta < 0) {
//                        delta = (sqrt(fabs(delta)) <= eps) ? 0 : sgn(c1)*sqrt(pow(c1,2)*pow(lam3,4));
                        delta = 0;
                    }
                    else {
                        delta = sgn(c1)*sqrt(delta);
                    }
                    f = (lam3*lam3*c1 + delta)/(2*pow(x,5));
                    x = lam3dp + dt*f;
                    double dx = dt*(f - fp);
                    pt.m_lam3d[i] = x;
                    if (fabs(dx) <= fabs(x)*errrel) cnvgd = true;
                    if (++iter == maxit) { error = true; }
                } while (!cnvgd && !error);
                if (error)
                    feLogWarning("SIV dashpot stretch calculation did not converge!");
            }
        }
    }
    
    for (int i=0; i<MAX_TERMS; ++i)
    {
        if ((m_t[i] > 0) && (m_g[i] > 0)) {
            double lams[3]; // Maxwell spring principal stretches
            double y = pt.m_lam3d[i];
            lams[0] = pt.m_lam[0]; lams[1] = pt.m_lam[1]; lams[2] = pt.m_lam[2]/y;
            mat3ds Us = pt.m_U[0]*lams[0] + pt.m_U[1]*lams[1] + pt.m_U[2]*lams[2];
            // evaluate Fs and Js to calculate stress in Maxwell spring
            ep.m_F = pt.m_R*Us; // Fs
            ep.m_J = Us.det();  // Js
            s += m_Mxwl->Stress(mp)*(m_g[i]*ep.m_J/Jsafe);
        }
    }
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;

    // return the total Cauchy stress,
    return s;
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FESIVViscoelastic::Tangent(FEMaterialPoint& mp)
{
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

    // store safe copy of deformation gradient
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;
    
    // get the viscoelastic point data
    FESIVViscoelasticMaterialPoint& pt = *mp.ExtractData<FESIVViscoelasticMaterialPoint>();

    // Calculate the new elastic Cauchy stress
    tens4ds c = m_Base->Tangent(mp)*m_g0;

    for (int i=0; i<MAX_TERMS; ++i)
    {
        if ((m_t[i] > 0) && (m_g[i] > 0)) {
            double lams[3]; // Maxwell spring principal stretches
            double y = pt.m_lam3d[i];
            lams[0] = pt.m_lam[0]; lams[1] = pt.m_lam[1]; lams[2] = pt.m_lam[2]/y;
            mat3ds Us = pt.m_U[0]*lams[0] + pt.m_U[1]*lams[1] + pt.m_U[2]*lams[2];
            // evaluate Fs and Js to calculate stress in Maxwell spring
            ep.m_F = pt.m_R*Us; // Fs
            ep.m_J = Us.det();  // Js
            c += m_Mxwl->Tangent(mp)*(m_g[i]*ep.m_J/Jsafe);
        }
    }
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;

    // return the total elastic tangent,
    return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density function
double FESIVViscoelastic::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    if (dt == 0) return 0;
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the viscoelastic point data
    FESIVViscoelasticMaterialPoint& pt = *mp.ExtractData<FESIVViscoelasticMaterialPoint>();

    // Calculate the new elastic Cauchy stress
    double sed = m_Base->StrainEnergyDensity(mp)*m_g0;

    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;

    // calculate new history variables
    // terms are accumulated in s
    for (int i=0; i<MAX_TERMS; ++i)
    {
        if ((m_t[i] > 0) && (m_g[i] > 0)) {
            double lams[3]; // Maxwell spring principal stretches
            double y = pt.m_lam3d[i];
            lams[0] = pt.m_lam[0]; lams[1] = pt.m_lam[1]; lams[2] = pt.m_lam[2]/y;
            mat3ds Us = pt.m_U[0]*lams[0] + pt.m_U[1]*lams[1] + pt.m_U[2]*lams[2];
            // evaluate Fs and Js to calculate stress in Maxwell spring
            ep.m_F = pt.m_R*Us; // Fs
            ep.m_J = Us.det();  // Js
            sed += m_Mxwl->StrainEnergyDensity(mp)*m_g[i];
        }
    }
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;

    // return the total Cauchy stress,
    return sed;
}
