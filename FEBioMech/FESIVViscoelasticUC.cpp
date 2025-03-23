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
#include "FESIVViscoelasticUC.h"
#include "FEUncoupledMaterial.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/DumpStream.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESIVViscoelasticUC, FEUncoupledMaterial)

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
    ADD_PARAMETER(m_ttype, FE_RANGE_CLOSED(0,1), "trigger");

    // define the material properties
    ADD_PROPERTY(m_Base, "parallel");
    ADD_PROPERTY(m_Mxwl, "Maxwell");

END_FECORE_CLASS();
//-----------------------------------------------------------------------------
//! constructor
FESIVViscoelasticUC::FESIVViscoelasticUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
    m_g0 = 1;
    for (int i=0; i<MAX_TERMS; ++i)
    {
        m_t[i] = 1;
        m_g[i] = 0;
    }
    m_ttype = 0;
    m_binit = false;

    m_Base = nullptr;
    m_Mxwl = nullptr;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPointData* FESIVViscoelasticUC::CreateMaterialPointData()
{
    return new FESIVViscoelasticMaterialPoint(m_Base->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! data initialization and checking
bool FESIVViscoelasticUC::Init()
{
    // combine bulk modulus from base material, Maxwell material and uncoupled viscoelastic material
    if (m_binit == false) m_K += m_Base->m_K + m_Mxwl->m_K;

    if (FEUncoupledMaterial::Init() == false) return false;

    m_binit = true;

    return true;
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FESIVViscoelasticUC::DevStress(FEMaterialPoint& mp)
{
    return mat3ds(0, 0, 0, 0, 0, 0);
/*    double dt = GetFEModel()->GetTime().timeIncrement;
    if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the viscoelastic point data
    FESIVViscoelasticMaterialPoint& pt = *mp.ExtractData<FESIVViscoelasticMaterialPoint>();

    // Calculate the new elastic Cauchy stress
    mat3ds s = m_Base->DevStress(mp)*m_g0;

    // evaluate right stretch, right Hencky, rotation
    mat3ds C = ep.RightCauchyGreen();
    double l2[3], l[3];
    vec3d v[3];
    C.eigen2(l2, v);
    l[0] = sqrt(l2[0]); l[1] = sqrt(l2[1]); l[2] = sqrt(l2[2]);
    mat3ds V[3];
    V[0] = dyad(v[0]); V[1] = dyad(v[1]); V[2] = dyad(v[2]);
    mat3ds U = V[0]*l[0] + V[1]*l[1] + V[2]*l[2];
    mat3ds Ui= V[0]/l[0] + V[1]/l[1] + V[2]/l[2];
    mat3ds H = V[0]*log(l[0]) + V[1]*log(l[1]) + V[2]*log(l[2]);
    mat3ds HydH = C.inverse()*((C*H).trace()/3);
    mat3ds DevH = H - HydH;
    mat3d R = ep.m_F*Ui;
    
    switch (m_ttype) {
        case 0:
        {
            // trigger in response to any strain
            pt.m_H = H;
        }
            break;
        case 1:
        {
            // trigger in response to distortional strain
            pt.m_H = DevH;
        }
            break;
            
        default:
            break;
    }
    mat3ds H_dot = (pt.m_H - pt.m_Hp)/dt;

    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;

    // calculate new history variables
    // terms are accumulated in s
    for (int i=0; i<MAX_TERMS; ++i)
    {
        if (m_t[i] > 0) {
            pt.m_Hs[i] = (pt.m_Hsp[i] + dt*H_dot)/(1+dt/m_t[i]);
            double den = pt.m_H.dotdot(pt.m_H);
            pt.m_alpha[i] = (den > 0) ? pt.m_Hs[i].dotdot(pt.m_H)/den : pt.m_alphap[i];
            if (pt.m_alpha[i] < 0) pt.m_alpha[i] = 0;
            if (pt.m_alpha[i] > 1) pt.m_alpha[i] = 1;
            mat3ds Ua = V[0]*pow(l[0],pt.m_alpha[i]) + V[1]*pow(l[1],pt.m_alpha[i]) + V[2]*pow(l[2],pt.m_alpha[i]);
            ep.m_F = R*Ua;
            ep.m_J = ep.m_F.det();
            s += m_Mxwl->DevStress(mp)*m_g[i];
        }
    }
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;

    // return the total Cauchy stress,
    return s;*/
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FESIVViscoelasticUC::DevTangent(FEMaterialPoint& mp)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    tens4ds c; c.zero();
/*    if (dt == 0) return c;
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the viscoelastic point data
    FESIVViscoelasticMaterialPoint& pt = *mp.ExtractData<FESIVViscoelasticMaterialPoint>();

    // Calculate the new elastic Cauchy stress
    c = m_Base->DevTangent(mp)*m_g0;

    // evaluate right stretch, right Hencky, rotation
    mat3ds C = ep.RightCauchyGreen();
    double l2[3], l[3];
    vec3d v[3];
    C.eigen2(l2, v);
    l[0] = sqrt(l2[0]); l[1] = sqrt(l2[1]); l[2] = sqrt(l2[2]);
    mat3ds V[3];
    V[0] = dyad(v[0]); V[1] = dyad(v[1]); V[2] = dyad(v[2]);
    mat3ds U = V[0]*l[0] + V[1]*l[1] + V[2]*l[2];
    mat3ds Ui= V[0]/l[0] + V[1]/l[1] + V[2]/l[2];
    mat3d R = ep.m_F*Ui;
    
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;

    // calculate new history variables
    // terms are accumulated in s
    for (int i=0; i<MAX_TERMS; ++i)
    {
        if (m_t[i] > 0) {
            mat3ds Ua = V[0]*pow(l[0],pt.m_alpha[i]) + V[1]*pow(l[1],pt.m_alpha[i]) + V[2]*pow(l[2],pt.m_alpha[i]);
            ep.m_F = R*Ua;
            ep.m_J = ep.m_F.det();
            c += m_Mxwl->DevTangent(mp)*m_g[i];
        }
    }
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;
*/
    // return the total Cauchy stress,
    return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density function
double FESIVViscoelasticUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    return 0;
/*    double dt = GetFEModel()->GetTime().timeIncrement;
    if (dt == 0) return 0;
    
    // get the elastic part
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the viscoelastic point data
    FESIVViscoelasticMaterialPoint& pt = *mp.ExtractData<FESIVViscoelasticMaterialPoint>();

    // Calculate the new elastic Cauchy stress
    double sed = m_Base->DevStrainEnergyDensity(mp)*m_g0;

    // evaluate right stretch, right Hencky, rotation
    mat3ds C = ep.RightCauchyGreen();
    double l2[3], l[3];
    vec3d v[3];
    C.eigen2(l2, v);
    l[0] = sqrt(l2[0]); l[1] = sqrt(l2[1]); l[2] = sqrt(l2[2]);
    mat3ds V[3];
    V[0] = dyad(v[0]); V[1] = dyad(v[1]); V[2] = dyad(v[2]);
    mat3ds U = V[0]*l[0] + V[1]*l[1] + V[2]*l[2];
    mat3ds Ui= V[0]/l[0] + V[1]/l[1] + V[2]/l[2];
    mat3d R = ep.m_F*Ui;
    
    mat3d Fsafe = ep.m_F;
    double Jsafe = ep.m_J;

    // calculate new history variables
    // terms are accumulated in s
    for (int i=0; i<MAX_TERMS; ++i)
    {
        if (m_t[i] > 0) {
            mat3ds Ua = V[0]*pow(l[0],pt.m_alpha[i]) + V[1]*pow(l[1],pt.m_alpha[i]) + V[2]*pow(l[2],pt.m_alpha[i]);
            ep.m_F = R*Ua;
            ep.m_J = ep.m_F.det();
            sed += m_Mxwl->DevStrainEnergyDensity(mp)*m_g[i];
        }
    }
    ep.m_F = Fsafe;
    ep.m_J = Jsafe;

    // return the total Cauchy stress,
    return sed;*/
}
