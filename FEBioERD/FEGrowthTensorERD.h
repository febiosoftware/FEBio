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
#include <FEBioMech/FEElasticFiberMaterial.h>
#include <FEBioMech/FEFiberMaterial.h>
#include <FEBioMech/febiomech_api.h>
#include "FEChemicalReactionERD.h"
#include "FEReactionERD.h"

//-----------------------------------------------------------------------------
//! Base class for growth tensors.
//!
class FEBIOERD_API FEGrowthTensorERD : public FEMaterialProperty
{
public:
    FECORE_BASE_CLASS(FEGrowthTensorERD);

public:
    FEGrowthTensorERD(FEModel* pfem);
    virtual ~FEGrowthTensorERD();

    //! growth tensor
    virtual mat3d GrowthTensor(FEMaterialPoint& pt, const vec3d& a0) = 0;

    //! inverse of growth tensor
    virtual mat3d GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& a0) = 0;

    //! referential solid density
    virtual double GrowthDensity(FEMaterialPoint& pt, const vec3d& a0) = 0;

    //! dFgdtheta
    virtual mat3ds dFgdtheta(FEMaterialPoint& pt, const vec3d& a0) = 0;

    //! dkdtheta
    virtual double dkdtheta(FEMaterialPoint& pt) = 0;

    //! dphidcdot
    virtual double dphidcdot(FEMaterialPoint& pt, double& sol_id) = 0;

    double SoluteConcentration(FEMaterialPoint& pt);

    double SBMConcentration(FEMaterialPoint& pt);

    vec3d UpdateNormal(FEMaterialPoint& pt, const vec3d& a0);

    double SpeciesGrowth(FEMaterialPoint& pt);

    double ActivationFunction(FEMaterialPoint& pt);

    double EnvironmentalFunction(FEMaterialPoint& pt);

    double Sigmoid(double a, double x, double x0, double b);

    double GrowthRate(FEMaterialPoint& pt);

    //! initialize
    bool Init() override;

public:
    FEParamVec3 m_fiber;
    FEParamDouble   m_gm;       //! isotropic growth multiplier
    //!SL: temporary place holder. Intent is to allow scaling by value of some state variable. In this case we are saying SBM id number to identify which SBM to base growth on. Will implement general solution later for more options.
    int             m_sbm_id = -1;   //! Which sbm should be used? Optional for now...
    int             m_sol_id = -1;
    double          theta_gamma = 0.2;
    double          theta_a = -2.0;
    double          k_min = 0.0;
    double          k_max = 1.0;
    //!SL: Flag to choose whether growth is based on current or referential concentration.
    bool m_referential_normal_flag = true;
protected:
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Volume growth
//!
class FEVolumeGrowthERD : public FEGrowthTensorERD
{
public:
    FEVolumeGrowthERD(FEModel* pfem) : FEGrowthTensorERD(pfem) { m_gm = 1; }
    virtual ~FEVolumeGrowthERD() {}

    //! growth tensor
    mat3d GrowthTensor(FEMaterialPoint& pt, const vec3d& a0) override;

    //! inverse of growth tensor
    mat3d GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& a0) override;

    //! referential solid density
    double GrowthDensity(FEMaterialPoint& pt, const vec3d& a0) override;

    //! dFgdtheta
    mat3ds dFgdtheta(FEMaterialPoint& pt, const vec3d& a0) override;

    //! dkdtheta
    virtual double dkdtheta(FEMaterialPoint& pt) override;

    //! dphidcdot
    virtual double dphidcdot(FEMaterialPoint& pt, double& sol_i) override;
};

//-----------------------------------------------------------------------------
//! Area growth
//!
class FEAreaGrowthERD : public FEGrowthTensorERD
{
public:
    FEAreaGrowthERD(FEModel* pfem) : FEGrowthTensorERD(pfem) { m_gm = 1; }
    virtual ~FEAreaGrowthERD() {}

    //! growth tensor
    mat3d GrowthTensor(FEMaterialPoint& pt, const vec3d& a0) override;

    //! inverse of growth tensor
    mat3d GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& a0) override;

    //! referential solid density
    double GrowthDensity(FEMaterialPoint& pt, const vec3d& a0) override;
    
    //! dFgdtheta
    mat3ds dFgdtheta(FEMaterialPoint& pt, const vec3d& a0) override;

    //! dkdtheta
    virtual double dkdtheta(FEMaterialPoint& pt) override;

    //! dphidcdot
    virtual double dphidcdot(FEMaterialPoint& pt, double& sol_i) override;

public:
    FEVec3dValuator* m_fiber_0;
    FEVec3dValuator* m_fiber_1;
};

//-----------------------------------------------------------------------------
//! Fiber growth
//!
class FEFiberGrowthERD : public FEGrowthTensorERD
{
public:
    FEFiberGrowthERD(FEModel* pfem) : FEGrowthTensorERD(pfem) { m_gm = 1; }
    virtual ~FEFiberGrowthERD() {}

    //! growth tensor
    mat3d GrowthTensor(FEMaterialPoint& pt, const vec3d& a0) override;

    //! inverse of growth tensor
    mat3d GrowthTensorInverse(FEMaterialPoint& pt, const vec3d& a0) override;

    //! referential solid density
    double GrowthDensity(FEMaterialPoint& pt, const vec3d& a0) override;

    //! dFgdtheta
    mat3ds dFgdtheta(FEMaterialPoint& pt, const vec3d& a0) override;

    //! dkdtheta
    virtual double dkdtheta(FEMaterialPoint& pt) override;

    //! dphidcdot
    virtual double dphidcdot(FEMaterialPoint& pt, double& sol_i) override;
};