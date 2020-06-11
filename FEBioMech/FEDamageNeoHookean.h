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
#include "FEDamageMaterialPoint.h"

//-----------------------------------------------------------------------------
// This material is a first attempt to include damage in hyper-elastic materials.
// It assumes the simple damage model as defined in Simo, CMAME 60 (1987), 153-173

//-----------------------------------------------------------------------------
class FEDamageNeoHookean : public FEElasticMaterial
{
public:
	FEDamageNeoHookean(FEModel* pfem);

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

	double	m_alpha;	//!< damage parameter alpha
	double	m_beta;		//!< damage parameter beta

protected:
	double	m_lam;
	double	m_mu;

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	//! data initialization and checking
	bool Init() override;

	// returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() override
	{ 
		return new FEDamageMaterialPoint(new FEElasticMaterialPoint);
	}

	// calculate damage reduction factor
	double Damage(FEMaterialPoint& pt);

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
