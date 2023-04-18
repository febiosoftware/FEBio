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

//-----------------------------------------------------------------------------
//! Growth tensor
//!
// define the material parameters
BEGIN_FECORE_CLASS(FEGrowthTensor, FEMaterialProperty)
    ADD_PROPERTY(m_fiber, "fiber", FEProperty::Optional)->SetDefaultType("vector");
END_FECORE_CLASS();

FEGrowthTensor::FEGrowthTensor(FEModel* pfem) : FEMaterialProperty(pfem) 
{ 
	m_fiber = nullptr; 
}

FEGrowthTensor::~FEGrowthTensor() {}

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

// define the material parameters
BEGIN_FECORE_CLASS(FEVolumeGrowth, FEGrowthTensor)
    ADD_PARAMETER(m_gm, "multiplier")->setLongName("volumetric growth multiplier");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! growth tensor
mat3d FEVolumeGrowth::GrowthTensor(FEMaterialPoint& pt, const vec3d& n0)
{
    return mat3dd(m_gm(pt));
}

//-----------------------------------------------------------------------------
//! inverse of growth tensor
mat3d FEVolumeGrowth::GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& n0)
{
    return mat3dd(1./m_gm(pt));
}

//-----------------------------------------------------------------------------
//! referential solid density
double FEVolumeGrowth::GrowthDensity(FEMaterialPoint& pt, const vec3d& n0)
{
    return pow(m_gm(pt), 3);
}

//-----------------------------------------------------------------------------
//! Area growth
//!

// define the material parameters
BEGIN_FECORE_CLASS(FEAreaGrowth, FEGrowthTensor)
    ADD_PARAMETER(m_gm , "multiplier")->setLongName("areal growth multiplier");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! growth tensor
mat3d FEAreaGrowth::GrowthTensor(FEMaterialPoint& pt, const vec3d& n0)
{
    double gmiso = sqrt(m_gm(pt));
    double gmani = 2 - gmiso;
    mat3d Fg = mat3dd(gmiso) + (n0 & n0)*(gmani - 1);
    return Fg;
}

//-----------------------------------------------------------------------------
//! inverse of growth tensor
mat3d FEAreaGrowth::GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& n0)
{
    double gmiso = sqrt(m_gm(pt));
    double gmani = 2 - gmiso;
    mat3d Fgi = mat3dd(1./gmiso) - (n0 & n0)*((gmani - 1)/gmiso/(gmiso+gmani-1));
    return Fgi;
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

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberGrowth, FEGrowthTensor)
    ADD_PARAMETER(m_gm , "multiplier")->setLongName("uniaxial growth multiplier");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! growth tensor
mat3d FEFiberGrowth::GrowthTensor(FEMaterialPoint& pt, const vec3d& n0)
{
    double gmiso = 1;
    double gmani = m_gm(pt);
    mat3d Fg = mat3dd(gmiso) + (n0 & n0)*(gmani - 1);
    return Fg;
}

//-----------------------------------------------------------------------------
//! inverse of growth tensor
mat3d FEFiberGrowth::GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& n0)
{
    double gmiso = 1;
    double gmani = m_gm(pt);
    mat3d Fgi = mat3dd(1./gmiso) - (n0 & n0)*((gmani - 1)/gmiso/(gmiso+gmani-1));
    return Fgi;
}

//-----------------------------------------------------------------------------
//! referential solid density
double FEFiberGrowth::GrowthDensity(FEMaterialPoint& pt, const vec3d& n0)
{
    return m_gm(pt);
}
