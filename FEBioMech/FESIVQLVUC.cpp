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
    mat3ds s = m_Base->Stress(mp);
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the viscoelastic point data
    FESIVQLVMaterialPoint& pt = *mp.ExtractData<FESIVQLVMaterialPoint>();
    
    // evaluate total strain
    mat3dd I(1);
    mat3ds C = ep.DevRightCauchyGreen();
    double lam2[3], lam[3];
    vec3d u[3];
    mat3ds U[3];
    double Ev[3];
    bool skip = false;
    // check if there is no (total) deformation
    if ((C-I).norm() <= eps) skip = true;
    if (!skip) {
        C.eigen2(lam2, u);
        // sort by strain magnitude
        Ev[0] = (lam2[0]-1)/2; Ev[1] = (lam2[1]-1)/2; Ev[2] = (lam2[2]-1)/2;
        lam[0] = sqrt(lam2[0]); lam[1] = sqrt(lam2[1]); lam[2] = sqrt(lam2[2]);
        U[0] = dyad(u[0]); U[1] = dyad(u[1]); U[2] = dyad(u[2]);
        // evaluate rotation tensor from inverse of stretch
        mat3ds Ui = U[0]/lam[0] + U[1]/lam[1] + U[2]/lam[2];
        pt.m_R = ep.m_F*Ui;
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
            //        reorder = true;
        }
        else return mat3ds(0, 0, 0, 0, 0, 0);
        
        // check and restore right-handedness of eigenvectors of C (if needed)
        double rhcs = (pt.m_u[0] ^ pt.m_u[1])*pt.m_u[2];
        if (rhcs < 0) {
            vec3d utmp = pt.m_u[0];
            double ltmp = pt.m_lam[0];
            mat3ds Utmp = pt.m_U[0];
            pt.m_u[0] = pt.m_u[1];
            pt.m_lam[0] = pt.m_lam[1];
            pt.m_U[0] = pt.m_U[1];
            pt.m_u[1] = utmp;
            pt.m_lam[1] = ltmp;
            pt.m_U[1] = Utmp;
        }
        // check and restore flipping of eigenvectors of C (if needed)
        if (pt.m_u[2]*pt.m_up[2] < 0) {
            pt.m_u[0] = -pt.m_u[0];
            pt.m_u[1] = -pt.m_u[1];
            pt.m_u[2] = -pt.m_u[2];
        }
    }
    else {
        pt.m_lam[0] = pt.m_lam[1] = pt.m_lam[2] = 1;
        // don't update pt.m_u or pt.m_U
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
    // force u3dot to be orthogonal with u3
    double e = 1 - pt.m_u[2]*pt.m_up[2];
    vec3d du3 = pt.m_u[2]*(1-e) - pt.m_up[2];
    vec3d u3dot = (du3.norm() > 1e-3) ? du3/dt : vec3d(0,0,0);
    
    double u3dot2 = u3dot*u3dot;
    mat3ds U3dot = dyads(u3dot, pt.m_u[2]);
    double lams[3]; // Maxwell spring principal stretches
    bool cnvgd = false;
    bool error = false;
    double f = 0;
    double fp = 0;
    double lam3 = pt.m_lam[2];
    double lam3dp = pt.m_lam3dp;
    do {
        fp = f;
        double x = pt.m_lam3d;
        double rx = sqrt(x);
        double x2 = x*x;
        lams[0] = pt.m_lam[0]*rx; lams[1] = pt.m_lam[1]*rx; lams[2] = pt.m_lam[2]/x;
        mat3ds Us = pt.m_U[0]*lams[0] + pt.m_U[1]*lams[1] + pt.m_U[2]*lams[2];
        mat3ds Es = ((Us*Us).sym() - I)/2;
        ep.m_F = pt.m_R*Us;
        ep.m_J = Us.det();
        mat3ds Smhat = m_Mxwl->PK2Stress(mp,Es)/(3*m_eta);
        double b = Smhat.dotdot((pt.m_U[2]-I)/(2*lams[2])+pt.m_U[2]*pow(lams[2],2))*x;
        double c = (x2*x)/3*pow(x-1./x2,2)*u3dot2
        - x2/2*(1/rx*(pow(pt.m_lam[2],2)-1/pt.m_lam[2]) - (pow(lams[2],2)-1/lams[2]))*Smhat.dotdot(U3dot);
        double delta = pow(b,2) - 4*c;
        if (delta < 0) {
            // in principle delta should not be negative, but this can happen due to finite time increments
            delta = 0;
        }
        else {
            delta = sgn(b)*sqrt(delta);
        }
        f = (b + delta)/2;
        x = lam3dp + dt*f;
        double dx = dt*(f - fp);
        pt.m_lam3d = x;
        if (fabs(dx) <= fabs(x)*errrel) cnvgd = true;
        if (++iter == maxit) { error = true; }
    } while (!cnvgd && !error);
    if (error)
        feLogWarning("SIV UC dashpot stretch calculation did not converge!");

    double y = pt.m_lam3d;
    lams[0] = pt.m_lam[0]*sqrt(y); lams[1] = pt.m_lam[1]*sqrt(y); lams[2] = pt.m_lam[2]/y;
    mat3ds Us = pt.m_U[0]*lams[0] + pt.m_U[1]*lams[1] + pt.m_U[2]*lams[2];
    // evaluate Fs and Js to calculate stress in Maxwell spring
    ep.m_F = pt.m_R*Us; // Fs
    ep.m_J = Us.det();  // Js
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
    tens4ds c = m_Base->Tangent(mp);

    // Calculate the base Cauchy stress
    mat3ds s = m_Base->Stress(mp);
    
    double lams[3]; // Maxwell spring principal stretches
    double y = pt.m_lam3d;
    lams[0] = pt.m_lam[0]*sqrt(y); lams[1] = pt.m_lam[1]*sqrt(y); lams[2] = pt.m_lam[2]/y;
    mat3ds Us = pt.m_U[0]*lams[0] + pt.m_U[1]*lams[1] + pt.m_U[2]*lams[2];
    // evaluate Fs and Js to calculate stress in Maxwell spring
    ep.m_F = pt.m_R*Us; // Fs
    ep.m_J = Us.det();  // Js
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

    // calculate new history variables
    // terms are accumulated in s
    double lams[3]; // Maxwell spring principal stretches
    double y = pt.m_lam3d;
    lams[0] = pt.m_lam[0]*sqrt(y); lams[1] = pt.m_lam[1]*sqrt(y); lams[2] = pt.m_lam[2]/y;
    mat3ds Us = pt.m_U[0]*lams[0] + pt.m_U[1]*lams[1] + pt.m_U[2]*lams[2];
    // evaluate Fs and Js to calculate stress in Maxwell spring
    ep.m_F = pt.m_R*Us; // Fs
    ep.m_J = Us.det();  // Js
    sed += m_Mxwl->DevStrainEnergyDensity(mp);

    ep.m_F = Fsafe;
    ep.m_J = Jsafe;

    // return the total Cauchy stress,
    return sed;
}
