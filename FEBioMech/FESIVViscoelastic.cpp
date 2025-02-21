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
#include <limits>
#include <cmath>


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
    m_H = m_Hp = m_dHp = mat3ds(0);
    m_R = mat3dd(1);
    for (int i=0; i<MAX_TERMS; ++i) {
        m_alpha[i] = m_alphap[i] = 1.0;
        m_dalphap[i] = 0;
        m_Hs[i] = m_Hsp[i] = m_dHsp[i] = mat3ds(0);
        m_mumr[i] = 0;
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
    m_dHp = (m_H - m_Hp)/dt;
    m_Hp = m_H;
    
    // copy previous data
    for (int i=0; i<MAX_TERMS; ++i) {
        m_dHsp[i] = (m_Hs[i] - m_Hsp[i])/dt;
        m_dalphap[i] = (m_alpha[i] - m_alphap[i])/dt;
        m_Hsp[i] = m_Hs[i];
        m_alphap[i] = m_alpha[i];
    }
    
    // don't forget to call the base class
    FEMaterialPointData::Update(timeInfo);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FESIVViscoelasticMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
    ar & m_H & m_Hp;
    ar & m_sed & m_sedp;
    ar & m_alpha;
    ar & m_Hs & m_Hsp & m_dHsp;
    ar & m_R;
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
    double dt = GetFEModel()->GetTime().timeIncrement;
    // evaluate viscosities in reference configuration
    if ((tp.timeStep == 0) && (tp.currentIteration == 0)) {
        mat3dd I(1.);

        // get the viscosities in the reference configuration
        FESIVViscoelasticMaterialPoint& pt = *mp.ExtractData<FESIVViscoelasticMaterialPoint>();
        for (int i=0; i<MAX_TERMS; ++i) {
            if (m_g[i] > 0) {
                tens4ds Cr = m_Mxwl->Tangent(mp);
                double Kr = (Cr.dot(I)).dotdot(I)/9.;
                double mur = ((Cr.dot2(I)).dotdot(I)-(Cr.dot(I)).dotdot(I)/3.)/10.;
                pt.m_mumr[i] = mur;
            }
        }
        return;
    }
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FESIVViscoelastic::Stress(FEMaterialPoint& mp)
{
    FETimeInfo& tp = GetFEModel()->GetTime();
    double dt = tp.timeIncrement;
    if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
    bool first = false;
    if ((tp.timeStep == 0) && tp.currentIteration == 0) first = true;

    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the viscoelastic point data
    FESIVViscoelasticMaterialPoint& pt = *mp.ExtractData<FESIVViscoelasticMaterialPoint>();
    
    // evaluate right stretch, right Hencky, rotation
    mat3dd I(1);
    mat3ds C = ep.RightCauchyGreen();
    double l2[3], l[3];
    vec3d u[3];
    C.eigen2(l2, u);
    l[0] = sqrt(l2[0]); l[1] = sqrt(l2[1]); l[2] = sqrt(l2[2]);
    mat3ds Ub[3];
    Ub[0] = dyad(u[0]); Ub[1] = dyad(u[1]); Ub[2] = dyad(u[2]);
    mat3ds H = Ub[0]*log(l[0]) + Ub[1]*log(l[1]) + Ub[2]*log(l[2]);
    
    pt.m_H = H;
    mat3ds H_dot = (pt.m_H - pt.m_Hp)/dt;
    double eps = 12*std::numeric_limits<double>::epsilon();
//    double eps = 1e-12;
    double residual = std::fmod(GetFEModel()->GetTime().currentTime,0.5);
    double H2 = H.dotdot(H);
    double HdH = H_dot.dotdot(H);
    double Hd2 = H_dot.dotdot(H_dot);
    double Hn = pt.m_H.norm();
    
    // Calculate the base PK2 stress
    mat3ds E = ep.Strain();
    mat3ds S = m_Base->PK2Stress(mp,E)*m_g0;
    mat3ds U = Ub[0]*l[0] + Ub[1]*l[1] + Ub[2]*l[2];
    mat3ds Ui= Ub[0]/l[0] + Ub[1]/l[1] + Ub[2]/l[2];
    mat3d R = ep.m_F*Ui;
    pt.m_R = R;
    
    // store safe copy of deformation gradient
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;
    
    // calculate new history variables
    // terms are accumulated in s
    double errrel = 1e-6;
    double errabs = 1e-12;
    int maxit = 100;
    int iter = 0;
    for (int i=0; i<MAX_TERMS; ++i)
    {
        // only solve for alpha if there is a non-zero time constant and gamma associated with this term
        if ((m_t[i] > 0) && (m_g[i] > 0)) {
            // if the strain is non-negligible, solve for alpha
//            if (((Hn > eps) && ((pt.m_alphap[i] >= -1./dt) || (pt.m_alphap[i] == 1./dt))) || first) {
//            if ((Hn > eps) || first) {
            if (Hn > eps) {
                // solve for alpha iteratively, starting with initial guess
                double y = (pt.m_alphap[i] < -1.0/dt) ? 1./dt : pt.m_alpha[i];
                double Dy = 1e-3;
                bool cnvgd = false;
                do {
                    mat3ds Ua = Ub[0]*pow(l[0],y) + Ub[1]*pow(l[1],y) +Ub[2]*pow(l[2],y);
                    ep.m_F = R*Ua;
                    ep.m_J = ep.m_F.det();
                    mat3ds Shat = m_Mxwl->PK2Stress(mp, ep.Strain())/(2*pt.m_mumr[i]);
                    // evaluate dSy/dy
                    Ua = Ub[0]*pow(l[0],y+Dy) + Ub[1]*pow(l[1],y+Dy) +Ub[2]*pow(l[2],y+Dy);
                    ep.m_F = R*Ua;
                    ep.m_J = ep.m_F.det();
                    mat3ds dShat = (m_Mxwl->PK2Stress(mp, ep.Strain())/(2*pt.m_mumr[i]) - Shat)/Dy;
                    
                    double f = (H2 + HdH*dt)*y + Shat.dotdot(H)*(dt/m_t[i]) - H2*pt.m_alphap[i] - HdH*dt;
                    double df = H2 + HdH*dt + dShat.dotdot(H)*(dt/m_t[i]);
                    double dy = -f/df;
                    y += dy;
                    pt.m_alpha[i] = y;
//                    if ((fabs(dy) <= fabs(y)*errrel) || (fabs(f) <= errabs)) cnvgd = true;
                    if (fabs(dy) <= fabs(y)*errrel) cnvgd = true;
                    if (++iter == maxit) cnvgd = true;
                    if (iter == maxit) {
                        bool error = true;
                    }
                } while (!cnvgd);
                pt.m_Hs[i] = H*pt.m_alpha[i];
            }
//            else if (Hn <= eps) {
            else {
                pt.m_Hs[i] = pt.m_Hsp[i] + pt.m_dHsp[i]*dt;
                pt.m_alpha[i] = (first) ? 1.0 : pt.m_dHsp[i].dotdot(H_dot)/Hd2;
            }
        }
    }
    
    for (int i=0; i<MAX_TERMS; ++i)
    {
        if ((m_t[i] > 0) && (m_g[i] > 0)) {
            mat3ds Ua = mat3dd(1);
            // calculation when strain ≠ 0
            if (Hn > eps) {
                double a = pt.m_alpha[i];
                Ua = Ub[0]*pow(l[0],a) + Ub[1]*pow(l[1],a) + Ub[2]*pow(l[2],a);
            }
            // when strain = 0, calculation is based on extrapolation from previous time step
            else {
                Ua = H2U(pt.m_Hs[i]);
            }
            ep.m_F = R*Ua;
            ep.m_J = ep.m_F.det();
            S += m_Mxwl->PK2Stress(mp, ep.Strain())*m_g[i];
        }
    }
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;
    // push-forward Piola transformation to get Cauchy stress
    mat3ds s = (ep.m_F*S*ep.m_F.transpose()).sym()/ep.m_J;

    // return the total Cauchy stress,
    return s;
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FESIVViscoelastic::Tangent(FEMaterialPoint& mp)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    tens4dmm c; c.zero();
    if (dt == 0) return c.supersymm();

    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

    // store safe copy of deformation gradient
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;
    
    // get the viscoelastic point data
    FESIVViscoelasticMaterialPoint& pt = *mp.ExtractData<FESIVViscoelasticMaterialPoint>();

    // Calculate the new elastic Cauchy stress
    c = m_Base->MaterialTangent(mp, ep.Strain())*m_g0;

    for (int i=0; i<MAX_TERMS; ++i)
    {
        if ((m_t[i] > 0) && (m_g[i] > 0)) {
            mat3ds Ua = H2U(pt.m_Hs[i]);
            ep.m_F = pt.m_R*Ua;
            ep.m_J = ep.m_F.det();
            c += m_Mxwl->MaterialTangent(mp, ep.Strain())*m_g[i];
        }
    }
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;

    // return the total elastic tangent,
    return (c.pp(Fsafe).supersymm())/Jsafe;
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

    // evaluate right stretch, right Hencky, rotation
    mat3ds C = ep.RightCauchyGreen();
    double l2[3], l[3];
    vec3d u[3];
    C.eigen2(l2, u);
    l[0] = sqrt(l2[0]); l[1] = sqrt(l2[1]); l[2] = sqrt(l2[2]);
    mat3ds Ub[3];
    Ub[0] = dyad(u[0]); Ub[1] = dyad(u[1]); Ub[2] = dyad(u[2]);
    mat3ds U = Ub[0]*l[0] + Ub[1]*l[1] + Ub[2]*l[2];
    mat3ds Ui= Ub[0]/l[0] + Ub[1]/l[1] + Ub[2]/l[2];
    mat3d R = ep.m_F*Ui;
    
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;

    // calculate new history variables
    // terms are accumulated in s
    for (int i=0; i<MAX_TERMS; ++i)
    {
        if (m_t[i] > 0) {
            double a = pt.m_alpha[i];
            mat3ds Ua = Ub[0]*pow(l[0],a) + Ub[1]*pow(l[1],a) + Ub[2]*pow(l[2],a);
            ep.m_F = R*Ua;
            ep.m_J = ep.m_F.det();
            sed += m_Mxwl->StrainEnergyDensity(mp)*m_g[i];
        }
    }
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;

    // return the total Cauchy stress,
    return sed;
}

//-----------------------------------------------------------------------------
//! convert right Hencky strain H to right stretch tensor U
mat3ds FESIVViscoelastic::H2U(mat3ds& H)
{
    double l[3];
    vec3d u[3];
    H.eigen2(l,u);
    mat3ds U = dyad(u[0])*exp(l[0]) + dyad(u[1])*exp(l[1]) + dyad(u[2])*exp(l[2]);
    return U;
}

//-----------------------------------------------------------------------------
//! solve a system of equations based on a tensorial nonlinear equation
bool FESIVViscoelastic::SolvedY(tens4ds& dTdY, mat3ds& Y, mat3ds& dY)
{
    // convert 4th-order tensor to 6x6 matrix
    double D[6][6];
    tens4ds dTi = dTdY.inverse();
    dTi.extract(D);
    double rhs[6] = {-Y.xx(), -Y.yy(), -Y.zz(), -Y.xy(), -Y.yz(), -Y.xz()};
    double sln[6] = {0};
    for (int i=0; i<6; ++i)
        for (int j=0; j<6; ++j) sln[i] += D[i][j]*rhs[j];
    dY = mat3ds(sln[0], sln[1], sln[2], sln[3], sln[4], sln[5]);
    return true;
}
