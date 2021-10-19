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
#include <FEBioMech/FEElasticMaterial.h>
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! Material class for the spherical fiber distribution with
//! fiber modulus dependent on sbm referential density

class FEBIOMIX_API FESFDSBM : public FEElasticMaterial
{
public:
	FESFDSBM(FEModel* pfem) : FEElasticMaterial(pfem) { m_alpha = 0;}
	
	//! Initialization
	bool Init() override;
	
	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp) override;
	
	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp) override;
	
	// Strain energy density
	virtual double StrainEnergyDensity(FEMaterialPoint& mp) override;
	
	//! return fiber modulus
	double FiberModulus(double rhor) { return m_ksi0*pow(rhor/m_rho0, m_g);}
	
	// declare the parameter list
	DECLARE_FECORE_CLASS();
	
public:
	double	m_alpha;	// coefficient of exponential argument
	double	m_beta;		// power in power-law relation
	double	m_ksi0;		// ksi = ksi0*(rhor/rho0)^gamma
    double  m_rho0;     // rho0
    double  m_g;        // gamma
	int		m_sbm;      //!< global id of solid-bound molecule
	int		m_lsbm;     //!< local id of solid-bound molecule
	
	static int		m_nres;	// integration rule
	static double	m_cth[];
	static double	m_sth[];
	static double	m_cph[];
	static double	m_sph[];
	static double	m_w[];
};
