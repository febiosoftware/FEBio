/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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



#include "stdafx.h"
#include "FEPermHolmesMow.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEPermHolmesMow, FEHydraulicPermeability)
	ADD_PARAMETER(m_perm , FE_RANGE_GREATER_OR_EQUAL(0.0), "perm" );
	ADD_PARAMETER(m_M    , FE_RANGE_GREATER_OR_EQUAL(0.0), "M"    );
	ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPermHolmesMow::FEPermHolmesMow(FEModel* pfem) : FEHydraulicPermeability(pfem)
{
	m_perm = 1;
	m_M = m_alpha = 0;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
mat3ds FEPermHolmesMow::Permeability(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// relative volume
	double J = et.m_J;
	// referential solid volume fraction
	double phi0 = pt.m_phi0;
	
	// --- strain-dependent isotropic permeability ---
	
	return mat3dd(m_perm*pow((J-phi0)/(1.0-phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0));
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4dmm FEPermHolmesMow::Tangent_Permeability_Strain(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// relative volume
	double J = et.m_J;
	// referential solid volume fraction
	double phi0 = pt.m_phi0;
	
	mat3dd I(1);	// Identity
	
	double k0 = m_perm*pow((J-phi0)/(1.0-phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	double K0prime = (J*J*m_M+(J*(m_alpha+1)-phi0)/(J-phi0))*k0;
	mat3ds k0hat = I*K0prime;
	
	return dyad1mm(I,k0hat)-dyad4s(I)*(2*k0);
}
