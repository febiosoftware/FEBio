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

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterial : public FEElasticMaterial
{
public:
    FEElasticFiberMaterial(FEModel* pfem);

	FEMaterialPointData* CreateMaterialPointData() override;
    
	// get the fiber vector (in global coordinates)
	vec3d FiberVector(FEMaterialPoint& mp);

	// calculate stress in fiber direction a0
	virtual mat3ds FiberStress(FEMaterialPoint& mp, const vec3d& a0) = 0;

	// Spatial tangent
	virtual tens4ds FiberTangent(FEMaterialPoint& mp, const vec3d& a0) = 0;

	//! Strain energy density
	virtual double FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) = 0;

    // Set or clear pre-stretch, as needed in multigenerational materials (e.g., reactive viscoelasticity)
    void SetPreStretch(const mat3ds Us) { m_Us = Us; m_bUs = true; }
    void ResetPreStretch() { m_bUs = false; }
    vec3d FiberPreStretch(const vec3d a0);
    
private:
	// These are made private since fiber materials should implement the functions above instead. 
	// The functions can still be reached when a fiber material is used in an elastic mixture. 
	// In those cases the fiber vector is taken from the first column of Q. 
	mat3ds Stress(FEMaterialPoint& mp) final { return FiberStress(mp, FiberVector(mp)); }
	tens4ds Tangent(FEMaterialPoint& mp) final { return FiberTangent(mp, FiberVector(mp)); }
	double StrainEnergyDensity(FEMaterialPoint& mp) final { return FiberStrainEnergyDensity(mp, FiberVector(mp)); }

private:
    mat3ds  m_Us;   //!< pre-stretch tensor for fiber
    bool    m_bUs;  //!< flag for pre-stretch
    
public:
	FEVec3dValuator*	m_fiber;	//!< fiber orientation

	DECLARE_FECORE_CLASS();
};

// helper class for constructing elastic fiber materials from classes derived from FEFiberMaterial. 
template <class fiberMat> 
class FEElasticFiberMaterial_T : public FEElasticFiberMaterial
{
public:
	FEElasticFiberMaterial_T(FEModel* fem) : FEElasticFiberMaterial(fem), m_fib(fem) {}

	bool Init() override 
	{ 
		if (FEElasticFiberMaterial::Init() == false) return false;
		return m_fib.Init(); 
	}
	bool Validate() override 
	{ 
		if (FEElasticFiberMaterial::Validate() == false) return false;
		return m_fib.Validate();
	}
	FEMaterialPointData* CreateMaterialPointData() override 
	{ 
		return new FEElasticMaterialPoint(m_fib.CreateMaterialPointData());
	}
	void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override 
	{
		FEElasticFiberMaterial::UpdateSpecializedMaterialPoints(mp, tp);
		m_fib.UpdateSpecializedMaterialPoints(mp, tp);
	}

	mat3ds FiberStress(FEMaterialPoint& mp, const vec3d& a0) override { return m_fib.FiberStress(mp, a0); }
	tens4ds FiberTangent(FEMaterialPoint& mp, const vec3d& a0) override { return m_fib.FiberTangent(mp, a0); }
	double FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override { return m_fib.FiberStrainEnergyDensity(mp, a0); }

protected:
	fiberMat	m_fib;
};
