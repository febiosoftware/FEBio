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
#include "FEPermRefTransIso.h"
#include <FECore/log.h>


// define the material parameters
BEGIN_FECORE_CLASS(FEPermRefTransIso, FEHydraulicPermeability)
	ADD_PARAMETER(m_perm0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "perm0")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_perm1T, FE_RANGE_GREATER_OR_EQUAL(0.0), "perm1T")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_perm1A, FE_RANGE_GREATER_OR_EQUAL(0.0), "perm1A")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_perm2T, FE_RANGE_GREATER_OR_EQUAL(0.0), "perm2T")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_perm2A, FE_RANGE_GREATER_OR_EQUAL(0.0), "perm2A")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_M0    , FE_RANGE_GREATER_OR_EQUAL(0.0), "M0");
	ADD_PARAMETER(m_MT    , FE_RANGE_GREATER_OR_EQUAL(0.0), "MT");
	ADD_PARAMETER(m_MA    , FE_RANGE_GREATER_OR_EQUAL(0.0), "MA");
	ADD_PARAMETER(m_alpha0, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha0");
	ADD_PARAMETER(m_alphaT, FE_RANGE_GREATER_OR_EQUAL(0.0), "alphaT");
	ADD_PARAMETER(m_alphaA, FE_RANGE_GREATER_OR_EQUAL(0.0), "alphaA");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPermRefTransIso::FEPermRefTransIso(FEModel* pfem) : FEHydraulicPermeability(pfem)
{
	m_perm0 = 1;
	m_perm1T = m_perm1A = 0;
	m_perm2T = m_perm2A = 0;
	m_M0 = m_MT = m_MA = 0;
	m_alpha0 = m_alphaT = m_alphaA = 0;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
mat3ds FEPermRefTransIso::Permeability(FEMaterialPoint& mp)
{
	vec3d V;			// axial material directions in reference configuration
	mat3ds m;			// axial texture tensor in current configuration
	
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// deformation gradient
	mat3d &F = et.m_F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.m_J;
	// referential solid volume fraction
	double phi0 = pt.m_phi0;
    double phisr = pt.m_phi0t;
	
    // check for potential error
    if (J <= phisr) feLogError("The perm-ref-trans-iso permeability calculation failed!\nThe volume ratio (J=%g) dropped below its theoretical minimum phi0=%g.",J,phisr);
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// Copy the texture direction in the reference configuration to V
	V.x = Q[0][0]; V.y = Q[1][0]; V.z = Q[2][0];
	m = dyad(F*V);	// Evaluate texture tensor in the current configuration
	
	// --- strain-dependent permeability ---
	
	double f, k1T, k1A, k2T, k2A;
	double k0 = m_perm0*pow((J-phisr)/(1-phi0),m_alpha0)*exp(m_M0*(J*J-1.0)/2.0);
	// Transverse direction
	f = pow((J-phisr)/(1-phi0),m_alphaT)*exp(m_MT*(J*J-1.0)/2.0);
	k1T = m_perm1T/(J*J)*f;
	k2T = 0.5*m_perm2T/pow(J,4)*f;
	// Axial direction
	f = pow((J-phisr)/(1-phi0),m_alphaA)*exp(m_MA*(J*J-1.0)/2.0);
	k1A = m_perm1A/(J*J)*f;
	k2A = 0.5*m_perm2A/pow(J,4)*f;
	// Permeability
	mat3ds kt = k0*I + k1T*b + (k1A-k1T)*m + 2*k2T*b.sqr() + (m*b).sym()*(2.0*(k2A-k2T));
	
	return kt;
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4dmm FEPermRefTransIso::Tangent_Permeability_Strain(FEMaterialPoint &mp)
{
	vec3d V;			// axial material directions in reference configuration
	mat3ds m;			// axial texture tensor in current configuration
	
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// Identity
	mat3dd I(1);
	
	// deformation gradient
	mat3d &F = et.m_F;
	
	// left cauchy-green matrix
	mat3ds b = et.LeftCauchyGreen();
	
	// relative volume
	double J = et.m_J;
	// referential solid volume fraction
	double phi0 = pt.m_phi0;
    double phisr = pt.m_phi0t;
	
    // check for potential error
    if (J <= phisr) feLogError("The perm-ref-trans-iso permeability calculation failed!\nThe volume ratio (J=%g) dropped below its theoretical minimum phi0=%g.",J,phisr);
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// Copy the texture direction in the reference configuration to V
	V.x = Q[0][0]; V.y = Q[1][0]; V.z = Q[2][0];
	m = dyad(F*V);	// Evaluate texture tensor in the current configuration
	
	double f, k0, K0prime;
	mat3ds k0hat, k1hat, k2hat;
	k0 = m_perm0*pow((J-phisr)/(1-phi0),m_alpha0)*exp(m_M0*(J*J-1.0)/2.0);
	K0prime = (1+J*(m_alpha0/(J-phisr)+m_M0*J))*k0;
	k0hat = mat3dd(K0prime);
	tens4dmm K4 = dyad1mm(I,k0hat)-dyad4s(I)*(2*k0);
	// Transverse direction
	f = pow((J-phisr)/(1-phi0),m_alphaT)*exp(m_MT*(J*J-1.0)/2.0);
	double k1T = m_perm1T/(J*J)*f;
	double k2T = 0.5*m_perm2T/pow(J,4)*f;
	mat3ds k1hatT = mat3dd((J*J*m_MT+(J*(m_alphaT-1)+phi0)/(J-phisr))*k1T);
	mat3ds k2hatT = mat3dd((J*J*m_MT+(J*(m_alphaT-3)+3*phi0)/(J-phisr))*k2T);
	// Axial direction
	f = pow((J-phisr)/(1-phi0),m_alphaA)*exp(m_MA*(J*J-1.0)/2.0);
	double k1A = m_perm1A/(J*J)*f;
	double k2A = 0.5*m_perm2A/pow(J,4)*f;
	mat3ds k1hatA = mat3dd((J*J*m_MA+(J*(m_alphaA-1)+phi0)/(J-phisr))*k1A);
	mat3ds k2hatA = mat3dd((J*J*m_MA+(J*(m_alphaA-3)+3*phi0)/(J-phisr))*k2A);
	//  Tangent
	K4 += dyad1mm(b.sqr(),k2hatT)*2 + dyad4s(b)*(4*k2T) + dyad4s(m,b)*(2*(k2A-k2T))
	+ dyad1mm(b,k1hatT) + dyad1mm(m,k1hatA-k1hatT) + dyad1mm((m*b).sym()*2.0,k2hatA-k2hatT);
	
	return K4;
}
