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
// TODO: This was introduced so that the FEGenerationMaterial fits in the current
//       framework, which assumes that all features classes are derived from base-classes.
class FEBIOMECH_API FEGenerationBase : public FEMaterialProperty
{
public:
	FEGenerationBase(FEModel* fem);

	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);

	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt);

	// returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;

public:
	double	btime;	//!< generation birth time

	FEElasticMaterial* m_pMat;	//!< pointer to elastic material

	FECORE_BASE_CLASS(FEGenerationBase)
};

//-----------------------------------------------------------------------------
//! Material defining a single generation of a multi-generation material
class FEGenerationMaterial : public FEGenerationBase
{
public:
	FEGenerationMaterial(FEModel* pfem);
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// forward declaration of material class
class FEElasticMultigeneration;

//-----------------------------------------------------------------------------
//! Multigenerational material point.
//! First generation exists at t=0. Second, third, etc. generations appear at t>0.
//! This material point stores the inverse of the relative deformation gradient of
//! second, third, etc. generations.  These relate the reference configuration of 
//! each generation relative to the first generation.

class FEMultigenerationMaterialPoint : public FEMaterialPointArray
{
public:
    FEMultigenerationMaterialPoint();
		
	FEMaterialPointData* Copy();

	//! data serialization
	void Serialize(DumpStream& ar);

	void Init();

	void Update(const FETimeInfo& timeInfo);

public:
	// multigenerational material data
	double	m_tgen;		//!< last generation time
    int     m_ngen;     //!< number of active generations
	FEElasticMultigeneration*	m_pmat;
};

//-----------------------------------------------------------------------------
//! Multigenerational solid

class FEElasticMultigeneration : public FEElasticMaterial
{
public:
	FEElasticMultigeneration(FEModel* pfem);
		
	// returns a pointer to a new material point object
    FEMaterialPointData* CreateMaterialPointData() override;

    // return number of materials
    int Materials() { return (int)m_MG.size(); }
    
    // return a generation material component
	FEGenerationBase* GetMaterial(int i) { return m_MG[i]; }
    
public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
		
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;
		
	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	int CheckGeneration(const double t);

public:
    double StrongBondSED(FEMaterialPoint& pt) override;
    double WeakBondSED(FEMaterialPoint& pt) override;
    
public:
	std::vector<FEGenerationBase*>	m_MG;		//!< multigeneration data

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
