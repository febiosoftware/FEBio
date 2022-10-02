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
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterialUC : public FEUncoupledMaterial
{
public:
    FEElasticFiberMaterialUC(FEModel* pfem);

	// Get the fiber direction (in global coordinates) at a material point
	vec3d FiberVector(FEMaterialPoint& mp);

	// calculate stress in fiber direction a0
	virtual mat3ds DevFiberStress(FEMaterialPoint& mp, const vec3d& a0) = 0;

	// Spatial tangent
	virtual tens4ds DevFiberTangent(FEMaterialPoint& mp, const vec3d& a0) = 0;

	//! Strain energy density
	virtual double DevFiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) = 0;
    
    // Set or clear pre-stretch, as needed in multigenerational materials (e.g., reactive viscoelasticity)
    void SetPreStretch(const mat3ds& Us) { m_Us = Us; m_bUs = true; }
    void ResetPreStretch() { m_bUs = false; }
    vec3d FiberPreStretch(const vec3d& a0);

private:
	// These are made private since fiber materials should implement the functions above instead. 
	// The functions can still be reached when a fiber material is used in an elastic mixture. 
	// In those cases the fiber vector is taken from the first column of Q. 
	mat3ds DevStress(FEMaterialPoint& mp) final { return DevFiberStress(mp, FiberVector(mp)); }
	tens4ds DevTangent(FEMaterialPoint& mp) final { return DevFiberTangent(mp, FiberVector(mp)); }
	double DevStrainEnergyDensity(FEMaterialPoint& mp) final { return DevFiberStrainEnergyDensity(mp, FiberVector(mp)); }
    
private:
    mat3ds  m_Us;   //!< pre-stretch tensor for fiber
    bool    m_bUs;  //!< flag for pre-stretch

public:
	
	FEVec3dValuator*	m_fiber;	//!< fiber orientation

	double	m_epsf;

	DECLARE_FECORE_CLASS();
};


template <class FiberMatUC>
class FEElasticFiberMaterialUC_T : public FEElasticFiberMaterialUC
{
public: 
	FEElasticFiberMaterialUC_T(FEModel* fem) : FEElasticFiberMaterialUC(fem), m_fib(fem) {}

	bool Init() override { return m_fib.Init(); }
	bool Validate() override { return m_fib.Validate(); }
	FEMaterialPointData* CreateMaterialPointData() override
	{
		FEMaterialPointData* mp = FEUncoupledMaterial::CreateMaterialPointData();
		mp->SetNext(m_fib.CreateMaterialPointData());
		return mp;
	}	
	void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override 
	{ 
		FEElasticFiberMaterialUC::UpdateSpecializedMaterialPoints(mp, tp);
		m_fib.UpdateSpecializedMaterialPoints(mp, tp); 
	}

	mat3ds DevFiberStress(FEMaterialPoint& mp, const vec3d& a0) override { return m_fib.DevFiberStress(mp, a0); }
	tens4ds DevFiberTangent(FEMaterialPoint& mp, const vec3d& a0) override { return m_fib.DevFiberTangent(mp, a0); }
	double DevFiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override { return m_fib.DevFiberStrainEnergyDensity(mp, a0); }

protected:
	FiberMatUC	m_fib;
};
