/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include <FEBioMix/febiomix_api.h>
#include "FEGrowthTensor.h"
#include <FECore/FEConstValueVec3.h>
#include <FECore/FEModel.h>
#include <FEBioMix/FEMultiphasic.h>
#include <FEBioMix/FEElasticReactionDiffusion.h>
#include <FEBioMix/FEChemicalReactionERD.h>
#include <FEBioMix/FEReactionERD.h>
#include <FECore/log.h>
#include <iostream>
#include <limits>

//-----------------------------------------------------------------------------
//! Growth tensor
//!
// define the material parameters
BEGIN_FECORE_CLASS(FEGrowthTensor, FEMaterialProperty)
ADD_PARAMETER(m_fiber, "fiber");
ADD_PARAMETER(m_gm, "multiplier")->setLongName("growth_multiplier");
ADD_PARAMETER(m_sbm_id, "sbm_id")->setLongName("sbm id for scaling growth");
ADD_PARAMETER(m_sol_id, "sol_id")->setLongName("sol id for scaling growth");
ADD_PARAMETER(m_referential_normal_flag, "referential_normal_direction");
ADD_PARAMETER(theta_gamma, "theta_gamma");
ADD_PARAMETER(theta_a, "theta_a");
ADD_PARAMETER(k_min, "k_min");
ADD_PARAMETER(k_max, "k_max");
END_FECORE_CLASS();

FEGrowthTensor::FEGrowthTensor(FEModel* pfem) : FEMaterialProperty(pfem)
{
    m_fiber = vec3d(1, 0, 0);
}

FEGrowthTensor::~FEGrowthTensor() {}

bool FEGrowthTensor::Init()
{
    return FEMaterialProperty::Init();
}

//-----------------------------------------------------------------------------
//! Volume growth
//!

//-----------------------------------------------------------------------------
//! growth tensor
mat3d FEVolumeGrowth::GrowthTensor(FEMaterialPoint& pt, const vec3d& n0)
{
    FEElasticMaterialPoint* ep = pt.ExtractData<FEElasticMaterialPoint>();
    FEKinematicMaterialPoint* kp = pt.ExtractData<FEKinematicMaterialPoint>();
    double theta = kp->m_theta;
    double gmiso = 1.0 + theta;
    return mat3dd(gmiso);
}

//-----------------------------------------------------------------------------
//! inverse of growth tensor
mat3d FEVolumeGrowth::GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& n0)
{
    mat3d Fg = GrowthTensor(pt, n0);
    return Fg.inverse();
}

//-----------------------------------------------------------------------------
//! referential solid density
double FEVolumeGrowth::GrowthDensity(FEMaterialPoint& pt, const vec3d& n0)
{
    /*mat3d Fg = GrowthTensor(pt, n0);
    return pow(Fg[0][0], 3);*/
    FEKinematicMaterialPoint& kp = *pt.ExtractData<FEKinematicMaterialPoint>();
    return pow(kp.m_theta, 3);
}

//-----------------------------------------------------------------------------
//! returns dFgdtheta for volume-type growth
mat3ds FEVolumeGrowth::dFgdtheta(FEMaterialPoint& pt, const vec3d& n0)
{
    double theta = this->GrowthRate(pt);
    return (1.0 / 3.0) * pow(theta, -2.0 / 3.0) * mat3dd(1.0);
}

//-----------------------------------------------------------------------------
//! returns dkdtheta for volume-type growth
double FEVolumeGrowth::dkdtheta(FEMaterialPoint& pt)
{
    //double theta = GrowthRate(pt);
    //double lhnum = exp((-theta_a - theta) / theta_gamma);
    //double rhnum = exp((theta_a - theta) / theta_gamma);
    //double lhden = theta_gamma * pow((lhnum + 1.0), 2.0);
    //double rhden = theta_gamma * pow((rhnum + 1.0), 2.0);

    //double dkdtheta = k_max * ((lhnum / lhden) + (rhnum / rhden));
    //return dkdtheta;
    return 0.0;
}

//-----------------------------------------------------------------------------
//! Area growth
//!

//-----------------------------------------------------------------------------
//! growth tensor
mat3d FEAreaGrowth::GrowthTensor(FEMaterialPoint& pt, const vec3d& n0)
{
    FEElasticMaterialPoint* ep = pt.ExtractData<FEElasticMaterialPoint>();
    FEKinematicMaterialPoint* kp = pt.ExtractData<FEKinematicMaterialPoint>();
    double sqrttheta = pow(kp->m_theta,0.5);
    vec3d n = m_referential_normal_flag ? n0 : UpdateNormal(pt, n0);
    mat3d Fg = sqrttheta * mat3dd(1.0) + (n & n) * (1.0 - sqrttheta);
    return Fg;
}

//-----------------------------------------------------------------------------
//! inverse of growth tensor
mat3d FEAreaGrowth::GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& n0)
{
    mat3d Fg = GrowthTensor(pt, n0);
    return Fg.inverse();
}

//-----------------------------------------------------------------------------
//! referential solid density
double FEAreaGrowth::GrowthDensity(FEMaterialPoint& pt, const vec3d& n0)
{
    FEKinematicMaterialPoint& kp = *pt.ExtractData<FEKinematicMaterialPoint>();
    return kp.m_theta;
}

//-----------------------------------------------------------------------------
//! returns dFgdtheta for area-type growth
mat3ds FEAreaGrowth::dFgdtheta(FEMaterialPoint& pt, const vec3d& n0)
{
    FEElasticMaterialPoint* ep = pt.ExtractData<FEElasticMaterialPoint>();
    FEKinematicMaterialPoint* kp = pt.ExtractData<FEKinematicMaterialPoint>();
    double sqrttheta = pow(kp->m_theta,0.5);
    //return (mat3dd(1.0) - (n0 & n0)).sym();
    return (0.5 / sqrttheta) * (mat3dd(1.0) - (n0 & n0)).sym();
}

//-----------------------------------------------------------------------------
//! returns dkdtheta for area-type growth
double FEAreaGrowth::dkdtheta(FEMaterialPoint& pt)
{
    //double theta = GrowthRate(pt);
    //double lhnum = exp((-theta_a - theta) / theta_gamma);
    //double rhnum = exp((theta_a - theta) / theta_gamma);
    //double lhden = theta_gamma * pow((lhnum + 1.0), 2.0);
    //double rhden = theta_gamma * pow((rhnum + 1.0), 2.0);

    //double dkdtheta = k_max * ((lhnum / lhden) + (rhnum / rhden));
    //return dkdtheta;+
    return 0.0;
}

//-----------------------------------------------------------------------------
//! Fiber growth
//!

//-----------------------------------------------------------------------------
//! growth tensor
mat3d FEFiberGrowth::GrowthTensor(FEMaterialPoint& pt, const vec3d& n0)
{
    FEElasticMaterialPoint* ep = pt.ExtractData<FEElasticMaterialPoint>();
    FEKinematicMaterialPoint* kp = pt.ExtractData<FEKinematicMaterialPoint>();
    double theta = kp->m_theta;
    double gmiso = 1.0;
    double gmani = theta - 1.0;
    vec3d n = m_referential_normal_flag ? n0 : UpdateNormal(pt, n0);
    mat3d Fg = mat3dd(gmiso) + (n & n) * gmani;
    return Fg;
}

//-----------------------------------------------------------------------------
//! inverse of growth tensor
mat3d FEFiberGrowth::GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& n0)
{
    mat3d Fg = GrowthTensor(pt, n0);
    return Fg.inverse();
}

//-----------------------------------------------------------------------------
//! referential solid density
double FEFiberGrowth::GrowthDensity(FEMaterialPoint& pt, const vec3d& n0)
{
    //return m_gm(pt);
    FEKinematicMaterialPoint& kp = *pt.ExtractData<FEKinematicMaterialPoint>();
    return m_gm(pt) * kp.m_theta;
}

//-----------------------------------------------------------------------------
//! returns dFgdtheta for fiber-type growth
mat3ds FEFiberGrowth::dFgdtheta(FEMaterialPoint& pt, const vec3d& n0)
{
    return (n0 & n0).sym();
}

//-----------------------------------------------------------------------------
//! returns dkdtheta for fiber-type growth
double FEFiberGrowth::dkdtheta(FEMaterialPoint& pt)
{
    double theta = GrowthRate(pt);
    double lhnum = exp((-theta_a - theta) / theta_gamma);
    double rhnum = exp((theta_a - theta) / theta_gamma);
    double lhden = theta_gamma * pow((lhnum + 1.0), 2.0);
    double rhden = theta_gamma * pow((rhnum + 1.0), 2.0);

    double dkdtheta = k_max * ((lhnum / lhden) + (rhnum / rhden));
    return dkdtheta;
}

//-----------------------------------------------------------------------------
//! bandpass activation function k(theta)
double FEGrowthTensor::ActivationFunction(FEMaterialPoint& pt)
{
    FEKinematicMaterialPoint& kp = *(pt.ExtractData<FEKinematicMaterialPoint>());
    double theta = kp.m_theta;
    double s_act = Sigmoid(1.0, theta, theta_a, theta_gamma);
    double s_inh = Sigmoid(1.0, theta, -1.0 * theta_a, theta_gamma);
    double k_theta = k_min + k_max * (s_act - s_inh);
    return k_theta;
}

//-----------------------------------------------------------------------------
//! microenvironmental function phi(dC/dt)
//! 
double FEGrowthTensor::EnvironmentalFunction(FEMaterialPoint& pt)
{
    FEElement* m_el = pt.m_elem;
    int nint = m_el->GaussPoints();
    FEDomain* dom = dynamic_cast<FEDomain*>(m_el->GetMeshPartition());
    FESolutesMaterialPoint& sp = *(pt.ExtractData<FESolutesMaterialPoint>());
    double c_n = SpeciesGrowth(pt);
    double c_n_p = sp.m_crp[m_sol_id - 1];
    double dt = GetTimeInfo().timeIncrement;
    double phi = (c_n - c_n_p) / dt;
    return (dt > 0.0) ? phi : 0.0;
}

//-----------------------------------------------------------------------------
//! get solute concentration helper function
double FEGrowthTensor::SoluteConcentration(FEMaterialPoint& pt)
{
    FEElement* m_el = pt.m_elem;
    FEDomain& dom = dynamic_cast<FEDomain&>(*m_el->GetMeshPartition());
    FESoluteInterface* pm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
    if (pm)
        return pm->GetActualSoluteConcentration(pt, m_sol_id - 1);
    else {
        FEElasticReactionDiffusionInterface* pm = dynamic_cast<FEElasticReactionDiffusionInterface*>(dom.GetMaterial());
        return pm->GetActualSoluteConcentration(pt, m_sol_id - 1);
    }
}

//-----------------------------------------------------------------------------
//! get SBM concentration helper function
double FEGrowthTensor::SBMConcentration(FEMaterialPoint& pt)
{
    FEElement* m_el = pt.m_elem;
    int nint = m_el->GaussPoints();
    FEDomain* dom = dynamic_cast<FEDomain*>(m_el->GetMeshPartition());
    FESolutesMaterialPoint& pd = *(pt.ExtractData<FESolutesMaterialPoint>());
    FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom->GetMaterial());
    return pm->SBMConcentration(pt, m_sbm_id - 1);
}

double FEGrowthTensor::SpeciesGrowth(FEMaterialPoint& pt)
{
    FEElement* m_el = pt.m_elem;
    int nint = m_el->GaussPoints();
    // if the growth tensor doesn't depend on concentration...
    if (m_sbm_id < 0 && m_sol_id < 0)
    {
        return 1.0;
    }

    //If we have a sbm dependence
    else if (m_sbm_id > 0 && m_sol_id < 0)
    {
        return SBMConcentration(pt);
    }

    //If we have sol dependence
    else if (m_sol_id > 0 && m_sbm_id < 0)
    {
        return SoluteConcentration(pt);
    }

    //If things aren't set properly just return the identity for now.
    else
    {
        return 1.0;
    }
}

vec3d FEGrowthTensor::UpdateNormal(FEMaterialPoint& pt, const vec3d& n0)
{
    FEElasticMaterialPoint* ep = pt.ExtractData<FEElasticMaterialPoint>();
    vec3d n = ep->m_F * n0;
    return n.normalized();
}

double FEGrowthTensor::Sigmoid(double a, double x, double x0, double b)
{
    double s = a / (1.0 + exp(-(x - x0) / b));
    return s;
}

double FEGrowthTensor::GrowthRate(FEMaterialPoint& pt)
{
    FEKinematicMaterialPoint* kp = pt.ExtractData<FEKinematicMaterialPoint>();
    double k_theta = ActivationFunction(pt);
    double phi = EnvironmentalFunction(pt);
    double dtheta = m_gm(pt) * k_theta * phi;
    double dt = this->CurrentTimeIncrement();    

    return dtheta;
}