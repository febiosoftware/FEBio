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

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

//-----------------------------------------------------------------------------
//! Material class for the ellipsoidal fiber distribution
//!

class FEEllipsoidalFiberDistribution : public FEElasticMaterial
{
	enum { MAX_INT = 45 };

public:
	FEEllipsoidalFiberDistribution(FEModel* pfem) : FEElasticMaterial(pfem) {}

	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp) override;

	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp) override;
	
	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
public:
	FEParamDouble	m_beta[3];	// power in power-law relation
	FEParamDouble	m_ksi[3];	// coefficient in power-law relation

								// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Material class for the ellipsoidal fiber distribution
//! (this is the old, obsolete implementation)

class FEEllipsoidalFiberDistributionOld : public FEElasticMaterial
{
public:
	FEEllipsoidalFiberDistributionOld(FEModel* pfem);
	
	//! Initialization
	bool Init() override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp) override;

	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp) override;
	
	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;

protected:
	void InitIntegrationRule();
    
public: // parameters
	FEParamDouble	m_beta[3];	// power in power-law relation
	FEParamDouble	m_ksi[3];	// coefficient in power-law relation

private:
	int		m_nres;	// integration rule
	double	m_cth[NSTH];
	double	m_sth[NSTH];
	double	m_cph[NSTH];
	double	m_sph[NSTH];
	double	m_w[NSTH];

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
