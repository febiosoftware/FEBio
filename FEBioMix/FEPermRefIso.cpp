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



#include "stdafx.h"
#include <FECore/log.h>
#include "FEPermRefIso.h"
#include <FECore/log.h>

#define INRANGE(x, a, b) ((x)>=(a) && (x)<=(b))

// define the material parameters
BEGIN_FECORE_CLASS(FEPermRefIso, FEHydraulicPermeability)
	ADD_PARAMETER(m_perm0, FE_RANGE_GREATER_OR_EQUAL(0.0), "perm0")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_perm1, FE_RANGE_GREATER_OR_EQUAL(0.0), "perm1")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_perm2, FE_RANGE_GREATER_OR_EQUAL(0.0), "perm2")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_M    , FE_RANGE_GREATER_OR_EQUAL(0.0), "M"    );
	ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPermRefIso::FEPermRefIso(FEModel* pfem) : FEHydraulicPermeability(pfem)
{
	m_perm0 = 1;
	m_perm1 = 0;
	m_perm2 = 0;
	m_M = m_alpha = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
bool FEPermRefIso::Validate()
{
	if (FEHydraulicPermeability::Validate() == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
mat3ds FEPermRefIso::Permeability(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.m_J;
	// referential solid volume fraction
	double phisr = pt.m_phi0t;
	
    // check for potential error
    if (J <= phisr) feLogError("The perm-ref-iso permeability calculation failed!\nThe volume ratio (J=%g) dropped below its theoretical minimum phi0=%g.",J,phisr);
    
	// --- strain-dependent permeability ---
	
	double f = pow((J-phisr)/(1-pt.m_phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	double k0 = m_perm0*f;
	double k1 = m_perm1/(J*J)*f;
	double k2 = 0.5*m_perm2/pow(J,4)*f;
	mat3ds kt = k0*I+k1*b+2*k2*b.sqr();
	
	return kt;
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4dmm FEPermRefIso::Tangent_Permeability_Strain(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.m_J;
	// referential solid volume fraction
	double phisr = pt.m_phi0t;
    double phi0 = pt.m_phi0;
	
    // check for potential error
    if (J <= phisr) feLogError("The perm-ref-iso permeability calculation failed!\nThe volume ratio (J=%g) dropped below its theoretical minimum phi0=%g.",J,phisr);
    
	double f = pow((J-phisr)/(1-phi0),m_alpha)*exp(m_M*(J*J-1.0)/2.0);
	double k0 = m_perm0*f;
	double k1 = m_perm1/(J*J)*f;
	double k2 = 0.5*m_perm2/pow(J,4)*f;
    double K0prime = (1+J*(m_alpha/(J-phi0)+m_M*J))*k0;
	double K1prime = (J*J*m_M+(J*(m_alpha-1)+phi0)/(J-phisr))*k1;
	double K2prime = (J*J*m_M+(J*(m_alpha-3)+3*phi0)/(J-phisr))*k2;
	mat3ds k0hat = I*K0prime;
	mat3ds k1hat = I*K1prime;
	mat3ds k2hat = I*K2prime;
	
	tens4dmm K4 = dyad1mm(I,k0hat)-dyad4s(I)*(2*k0)
	+ dyad1mm(b,k1hat)
	+ dyad1mm(b.sqr(),k2hat)*2+dyad4s(b)*(4*k2);
	
	return K4;
}
