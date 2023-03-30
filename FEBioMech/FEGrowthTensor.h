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

#pragma once
#include "FEElasticFiberMaterial.h"
#include "FEFiberMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for growth tensors.
//!
class FEGrowthTensor : public FEFiberMaterial
{
public:
    FEGrowthTensor(FEModel* pfem) : FEFiberMaterial(pfem) { m_fiber = nullptr; }
    virtual ~FEGrowthTensor(){}
    
    //! growth tensor
    virtual mat3d GrowthTensor(FEMaterialPoint& pt, const vec3d& a0) = 0;
    
    //! inverse of growth tensor
    virtual mat3d GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& a0) = 0;
    
    //! referential solid density
    virtual double GrowthDensity(FEMaterialPoint& pt, const vec3d& a0) = 0;
    
    //! stress
    mat3ds FiberStress(FEMaterialPoint& pt, const vec3d& a0) override { return mat3dd(0); }
    
    //! tangent
    tens4ds FiberTangent(FEMaterialPoint& pt, const vec3d& a0) override { tens4ds c; c.zero(); return c; }
    
    //! strain-energy density
    double FiberStrainEnergyDensity(FEMaterialPoint& pt, const vec3d& a0) override { return 0; }
    
    double SoluteConcentration(FEMaterialPoint& pt);

    double SBMConcentration(FEMaterialPoint& pt);

    double Multiplier(FEMaterialPoint& pt, const vec3d& a0) { };

    double SpeciesGrowth(FEMaterialPoint& pt);

    //! initialize
    bool Init() override;
    
public:
    FEVec3dValuator* m_fiber;
    FEParamDouble   m_gm;       //! isotropic growth multiplier
    //!SL: temporary place holder. Intent is to allow scaling by value of some state variable. In this case we are saying SBM id number to identify which SBM to base growth on. Will implement general solution later for more options.
    int             m_sbm_id = -1;   //! Which sbm should be used? Optional for now...
    int             m_sol_id = -1;

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Volume growth
//!
class FEVolumeGrowth : public FEGrowthTensor
{
public:
    FEVolumeGrowth(FEModel* pfem) : FEGrowthTensor(pfem) { m_gm = 1; }
    virtual ~FEVolumeGrowth(){}
    
    //! growth tensor
    mat3d GrowthTensor(FEMaterialPoint& pt, const vec3d& a0) override;
    
    //! inverse of growth tensor
    mat3d GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& a0) override;
    
    //! referential solid density
    double GrowthDensity(FEMaterialPoint& pt, const vec3d& a0) override;

public:
    
    // declare the parameter list
    //DECLARE_FECORE_CLASS();

};

class FEElasticVolumeGrowth : public FEElasticFiberMaterial_T<FEVolumeGrowth>
{
public:
    FEElasticVolumeGrowth(FEModel* fem) : FEElasticFiberMaterial_T<FEVolumeGrowth>(fem) {}
    //DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Area growth
//!
class FEAreaGrowth : public FEGrowthTensor
{
public:
    FEAreaGrowth(FEModel* pfem) : FEGrowthTensor(pfem) { m_gm = 1; }
    virtual ~FEAreaGrowth(){}
    
    //! growth tensor
    mat3d GrowthTensor(FEMaterialPoint& pt, const vec3d& a0) override;
    
    //! inverse of growth tensor
    mat3d GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& a0) override;
    
    //! referential solid density
    double GrowthDensity(FEMaterialPoint& pt, const vec3d& a0) override;
    
public:
    FEVec3dValuator* m_fiber_0;
    FEVec3dValuator* m_fiber_1;
    // declare the parameter list
    //DECLARE_FECORE_CLASS();

};

class FEElasticAreaGrowth : public FEElasticFiberMaterial_T<FEAreaGrowth>
{
public:
    FEElasticAreaGrowth(FEModel* fem) : FEElasticFiberMaterial_T<FEAreaGrowth>(fem) {}
    //DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Fiber growth
//!
class FEFiberGrowth : public FEGrowthTensor
{
public:
    FEFiberGrowth(FEModel* pfem) : FEGrowthTensor(pfem) { m_gm = 1; }
    virtual ~FEFiberGrowth(){}
    
    //! growth tensor
    mat3d GrowthTensor(FEMaterialPoint& pt, const vec3d& a0) override;
    
    //! inverse of growth tensor
    mat3d GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& a0) override;
    
    //! referential solid density
    double GrowthDensity(FEMaterialPoint& pt, const vec3d& a0) override;
    
public:

    // declare the parameter list
    // DECLARE_FECORE_CLASS();

};

class FEElasticFiberGrowth : public FEElasticFiberMaterial_T<FEFiberGrowth>
{
public:
    FEElasticFiberGrowth(FEModel* fem) : FEElasticFiberMaterial_T<FEFiberGrowth>(fem) {}
    //DECLARE_FECORE_CLASS();
};