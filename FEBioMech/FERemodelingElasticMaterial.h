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
//! Base class for solid supply.
//! These materials need to define the solid supply and tangent supply functions.
//! The solid supply has units of mass/(referential volume)/time
//!

class FEBIOMECH_API FESolidSupply : public FEMaterialProperty
{
public:
	//! constructor
	FESolidSupply(FEModel* pfem) : FEMaterialProperty(pfem) {}

	//! solid supply
	virtual double Supply(FEMaterialPoint& pt) = 0;
	
	//! tangent of solute supply with respect to strain
	virtual mat3ds Tangent_Supply_Strain(FEMaterialPoint& mp) = 0;
	
	//! tangent of solute supply with respect to referential density
	virtual double Tangent_Supply_Density(FEMaterialPoint& mp) = 0;	

	FECORE_BASE_CLASS(FESolidSupply)
};

//-----------------------------------------------------------------------------
//! Material point data for remodeling elastic materials
class FEBIOMECH_API FERemodelingMaterialPoint : public FEMaterialPointData
{
public:
	FERemodelingMaterialPoint(FEMaterialPointData*pt) : FEMaterialPointData(pt) {}
    
	FEMaterialPointData* Copy();
    
	void Init();
    
	void Update(const FETimeInfo& timeInfo);

	void Serialize(DumpStream& ar);

public:
	double		m_sed;		//!< strain energy density
	double		m_dsed;		//!< derivative of strain energy density with mass density
	double		m_rhor;		//!< current referential mass density
	double		m_rhorp;	//!< referential mass density at previous time step
};

//-----------------------------------------------------------------------------
//! A material that wants to use the remodeling framework needs to implement 
//! this additional interface.
class FERemodelingInterface
{
public:
	//! calculate strain energy density at material point
	virtual double StrainEnergy(FEMaterialPoint& pt) = 0;

	//! calculate tangent of strain energy density with solid density at material point
	virtual double Tangent_SE_Density(FEMaterialPoint& pt) = 0;

	//! calculate tangent of stress with solid density at material point
	virtual mat3ds Tangent_Stress_Density(FEMaterialPoint& pt) = 0;
};

//-----------------------------------------------------------------------------
//! Material class for remodeling solids
class FEBIOMECH_API FERemodelingElasticMaterial : public FEElasticMaterial
{
public:
	//! constructor
	FERemodelingElasticMaterial(FEModel* pfem);
	
	//! strain energy density function
	double StrainEnergyDensity(FEMaterialPoint& pt) override;
	
	//! stress function
	mat3ds Stress(FEMaterialPoint& pt) override;
	
	//! tangent function of stress with strain
	tens4ds Tangent(FEMaterialPoint& pt) override;
	
	//! tangent function of strain energy density with solid mass density
	double Tangent_SE_Density(FEMaterialPoint& pt);
	
	//! tangent function of stress with solid mass density
	mat3ds Tangent_Stress_Density(FEMaterialPoint& pt);
	
	// returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override
	{
		return new FERemodelingMaterialPoint(m_pBase->CreateMaterialPointData());
	}
    
	// get the elastic material
	FEElasticMaterial* GetElasticMaterial() override { return m_pBase; }

public:
	FEElasticMaterial*	m_pBase;		//!< pointer to elastic solid material
	FESolidSupply*		m_pSupp;		//!< pointer to solid supply material
	double				m_rhormin;		//!< minimum density
	double				m_rhormax;		//!< maximum density
	
public:
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
