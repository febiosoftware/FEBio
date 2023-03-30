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

#include "FEGrowthTensor.h"
#include <FECore/FEConstValueVec3.h>
#include <FECore/FEModel.h>
#include <FEBioMix/FEMultiphasic.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
//! Growth tensor
//!
// define the material parameters
BEGIN_FECORE_CLASS(FEGrowthTensor, FEFiberMaterial)
    ADD_PROPERTY(m_fiber, "fiber", FEProperty::Optional)->SetDefaultType("vector");
    ADD_PARAMETER(m_gm, "multiplier")->setLongName("time_multiplier");
    ADD_PARAMETER(m_sbm_id, "sbm_id")->setLongName("sbm id for scaling growth");
    ADD_PARAMETER(m_sol_id, "sol_id")->setLongName("sol id for scaling growth");
END_FECORE_CLASS();

bool FEGrowthTensor::Init()
{
    if (m_fiber == nullptr) {
        FEConstValueVec3* val = fecore_new<FEConstValueVec3>("vector", nullptr);
        val->value() = vec3d(1,0,0);
        m_fiber = val;
    }

    return true;
}

//-----------------------------------------------------------------------------
//! Volume growth
//!

//-----------------------------------------------------------------------------
//! growth tensor
mat3d FEVolumeGrowth::GrowthTensor(FEMaterialPoint& pt, const vec3d& n0)
{
    FEElement* m_e = pt.m_elem;
    FEDomain* dom = dynamic_cast<FEDomain*>(m_e->GetMeshPartition());
    // calculate the element volume
    FEMesh* mesh = dom->GetMesh();
    double g = 1.0 + m_gm(pt) * SpeciesGrowth(pt);
    return mat3dd(g);
}

//-----------------------------------------------------------------------------
//! inverse of growth tensor
mat3d FEVolumeGrowth::GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& n0)
{
    mat3d Fg = GrowthTensor(pt, n0);
    return Fg.inverse();
    //return mat3dd(1./m_gm(pt));
}

//-----------------------------------------------------------------------------
//! referential solid density
double FEVolumeGrowth::GrowthDensity(FEMaterialPoint& pt, const vec3d& n0)
{
    mat3d Fg = GrowthTensor(pt, n0);
    return pow(Fg[0][0], 3);
    //return pow(m_gm(pt), 3);
}

//-----------------------------------------------------------------------------
//! Area growth
//!

//-----------------------------------------------------------------------------
//! growth tensor
mat3d FEAreaGrowth::GrowthTensor(FEMaterialPoint& pt, const vec3d& n0)
{
    //! Gerard's implementation
    double gmiso = 1.0 + m_gm(pt) * SpeciesGrowth(pt);
    double gmani = 2.0 - gmiso;
    mat3d Fg = mat3dd(gmiso) + (n0 & n0) * (gmani - 1);
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
    return m_gm(pt);
}

//-----------------------------------------------------------------------------
//! Fiber growth
//!

//-----------------------------------------------------------------------------
//! growth tensor
mat3d FEFiberGrowth::GrowthTensor(FEMaterialPoint& pt, const vec3d& n0)
{
    double gmiso = 1;
    double gmani = 1.0 + m_gm(pt) * SpeciesGrowth(pt);
    mat3d Fg = mat3dd(gmiso) + (n0 & n0)*(gmani - 1);
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
    return m_gm(pt);
}

//-----------------------------------------------------------------------------
//! get solute concentration helper function
double FEGrowthTensor::SoluteConcentration(FEMaterialPoint& pt)
{
    FEElement* m_el = pt.m_elem;
    int nint = m_el->GaussPoints();
    FEDomain* dom = dynamic_cast<FEDomain*>(m_el->GetMeshPartition());
    FESolutesMaterialPoint& pd = *(pt.ExtractData<FESolutesMaterialPoint>());
    FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom->GetMaterial());
    return pd.m_c[m_sol_id - 1];
}

//-----------------------------------------------------------------------------
//! get SBM concentration helper function
double FEGrowthTensor::SBMConcentration(FEMaterialPoint & pt)
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
    // if the growth tensor doesn't depend on an sbm concentration...
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