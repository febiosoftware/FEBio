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
#pragma once
#include "FEElasticMaterial.h"
#include "FEElasticFiberMaterial.h"

// These classes implement the elastic damage framework from 
// Balzani, Brinkhues, Holzapfel, Comput. Methods Appl. Mech. Engrg. 213216 (2012) 139151

class FEFiberDamagePoint;

//=================================================================================================
// Base class for continuous damage elastic fiber materials. 
//
class FEDamageElasticFiber : public FEElasticFiberMaterial
{
public:
	FEDamageElasticFiber(FEModel* fem);

	double Damage(FEMaterialPoint& mp);
	double Damage(FEMaterialPoint& mp, int n);

	//! Strain energy density
	double FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;

	// calculate stress in fiber direction a0
	mat3ds FiberStress(FEMaterialPoint& mp, const vec3d& a0) override;

	// Spatial tangent
	tens4ds FiberTangent(FEMaterialPoint& mp, const vec3d& a0) override;

protected:
	// strain-energy and its derivatives
	virtual double Psi0(FEMaterialPoint& mp, const vec3d& a0);
	virtual mat3ds dPsi0_dC(FEMaterialPoint& mp, const vec3d& a0);
	virtual tens4ds d2Psi0_dC(FEMaterialPoint& mp, const vec3d& a0);

	virtual double m(double P);
	virtual double dm_dP(double P);
	virtual double d2m_dP(double P);

protected:
	// damage model parameters
	double	m_tinit;		// start time of damage
	double	m_Dmax;			// max damage
	double	m_beta_s;		// saturation parameter
	double	m_gamma_max;	// saturation parameter
	double	m_r_s, m_r_inf;

	// D2 parameters
	double m_D2_a;
	double m_D2_b;
	double m_D2_c;
	double m_D2_d;

	double m_D3_g0;
	double m_D3_rg;
	double m_D3_inf;

	DECLARE_FECORE_CLASS();
};

//=================================================================================================
class FEDamageFiberPower : public FEDamageElasticFiber
{
public:
	FEDamageFiberPower(FEModel* fem);

	FEMaterialPointData* CreateMaterialPointData() override;

protected:
	double Psi0(FEMaterialPoint& mp, const vec3d& a0) override;
	mat3ds dPsi0_dC(FEMaterialPoint& mp, const vec3d& a0) override;
	tens4ds d2Psi0_dC(FEMaterialPoint& mp, const vec3d& a0) override;

	double m(double P) override;
	double dm_dP(double P) override;
	double d2m_dP(double P) override;

public:
	// fiber parameters
	double	m_a1, m_a2, m_kappa;

	DECLARE_FECORE_CLASS();
};

//=================================================================================================
class FEDamageFiberExponential : public FEDamageElasticFiber
{
public:
	FEDamageFiberExponential(FEModel* fem);

	FEMaterialPointData* CreateMaterialPointData() override;

protected:
	double Psi0(FEMaterialPoint& mp, const vec3d& a0) override;
	mat3ds dPsi0_dC(FEMaterialPoint& mp, const vec3d& a0) override;
	tens4ds d2Psi0_dC(FEMaterialPoint& mp, const vec3d& a0) override;

	double m(double P) override;
	double dm_dP(double P) override;
	double d2m_dP(double P) override;

public:
	// fiber parameters
	double	m_k1, m_k2, m_kappa;

	DECLARE_FECORE_CLASS();
};


//=================================================================================================
class FEDamageFiberExpLinear: public FEDamageElasticFiber
{
public:
	FEDamageFiberExpLinear(FEModel* fem);

	FEMaterialPointData* CreateMaterialPointData() override;

protected:
	double Psi0(FEMaterialPoint& mp, const vec3d& a0) override;
	mat3ds dPsi0_dC(FEMaterialPoint& mp, const vec3d& a0) override;
	tens4ds d2Psi0_dC(FEMaterialPoint& mp, const vec3d& a0) override;

	double m(double P) override;
	double dm_dP(double P) override;
	double d2m_dP(double P) override;

public:
	double	m_c3;
	double	m_c4;
	double	m_c5;
	double	m_lamax;

	DECLARE_FECORE_CLASS();
};
