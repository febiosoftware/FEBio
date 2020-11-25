/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEElasticMaterial.h"
#include "FEElasticFiberMaterial.h"

// These classes implement the elastic damage framework from 
// Balzani, Brinkhues, Holzapfel, Comput. Methods Appl. Mech. Engrg. 213–216 (2012) 139–151

class FEFiberDamagePoint;

class FEDamageInterface
{
public:
	FEDamageInterface() {}
	virtual ~FEDamageInterface() {}
	double Damage(FEMaterialPoint& mp);
};

class FEDamageFiberPower : public FEElasticFiberMaterial, public FEDamageInterface
{
public:
	FEDamageFiberPower(FEModel* fem);

	FEMaterialPoint* CreateMaterialPointData() override;

	//! Strain energy density
	double FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;

	// calculate stress in fiber direction a0
	mat3ds FiberStress(FEMaterialPoint& mp, const vec3d& a0) override;

	// Spatial tangent
	tens4ds FiberTangent(FEMaterialPoint& mp, const vec3d& a0) override;

public:
	// fiber parameters
	double	m_a1, m_a2, m_kappa;

	// damage model parameters
	double	m_tinit;		// start time of damage
	double	m_Dmax;			// max damage
	double	m_beta_s;		// saturation parameter
	double	m_gamma_max;	// saturation parameter
	double	m_r_s, m_r_inf;

	DECLARE_FECORE_CLASS();
};

class FEDamageFiberExponential : public FEElasticFiberMaterial, public FEDamageInterface
{
public:
	FEDamageFiberExponential(FEModel* fem);

	FEMaterialPoint* CreateMaterialPointData() override;

	//! Strain energy density
	double FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;

	// calculate stress in fiber direction a0
	mat3ds FiberStress(FEMaterialPoint& mp, const vec3d& a0) override;

	// Spatial tangent
	tens4ds FiberTangent(FEMaterialPoint& mp, const vec3d& a0) override;

public:
	// fiber parameters
	double	m_k1, m_k2, m_kappa;

	// damage model parameters
	double	m_tinit;		// start time of damage
	double	m_Dmax;			// max damage
	double	m_beta_s;		// saturation parameter
	double	m_gamma_max;	// saturation parameter
	double	m_r_s, m_r_inf;

	DECLARE_FECORE_CLASS();
};

class FEBIOMECH_API FEContinuousElasticDamage : public FEElasticMaterial, public FEDamageInterface
{
public:
	FEContinuousElasticDamage(FEModel* fem);

	FEMaterialPoint* CreateMaterialPointData() override;

public:
	double StrainEnergyDensity(FEMaterialPoint& pt) override;

	mat3ds Stress(FEMaterialPoint& pt) override;

	tens4ds Tangent(FEMaterialPoint& mp) override;

private:
	mat3ds MatrixStress(FEMaterialPoint& mp);
	tens4ds MatrixTangent(FEMaterialPoint& mp);

private:
	// matrix parameters
	double		m_c1;
	double		m_k;
	FEParamVec3		m_fiber;

	// fiber parameters
	FEDamageFiberPower	m_fib;

	DECLARE_FECORE_CLASS();
};
