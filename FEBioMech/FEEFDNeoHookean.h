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
#include "FECore/FEMaterial.h"
#include "FENeoHookean.h"
#include "FEEllipsoidalFiberDistribution.h"

class FEEFDNeoHookean :	public FEElasticMaterial
{
public:
	FEEFDNeoHookean(FEModel* pfem);

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio
	double	m_ksi[3];	//!< ksi
	double	m_beta[3];	//!< beta

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	//! data initialization
	bool Init() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

public:

	FEEllipsoidalFiberDistribution	m_EFD;
	FENeoHookean					m_NH;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};


//----------------------------------------------------------------
// "old" version using full-sphere integration
class FEEFDNeoHookeanOld :	public FEElasticMaterial
{
public:
	FEEFDNeoHookeanOld(FEModel* pfem) : FEElasticMaterial(pfem), m_EFD(pfem), m_NH(pfem) {}

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	//! data initialization and checking
	bool Init() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

public:

	FEEllipsoidalFiberDistributionOld	m_EFD;
	FENeoHookean						m_NH;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
