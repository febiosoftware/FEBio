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
#include "FEPermRefOrtho.h"
#include <FECore/log.h>


// define the material parameters
BEGIN_FECORE_CLASS(FEPermRefOrtho, FEHydraulicPermeability)
	ADD_PARAMETER(m_perm0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "perm0")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_M0    , FE_RANGE_GREATER_OR_EQUAL(0.0), "M0");
	ADD_PARAMETER(m_alpha0, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha0");
	ADD_PARAMETER(m_perm1, 3, "perm1")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_perm2, 3, "perm2")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_M    , 3, "M");
	ADD_PARAMETER(m_alpha, 3, "alpha");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPermRefOrtho::FEPermRefOrtho(FEModel* pfem) : FEHydraulicPermeability(pfem)
{
	m_perm0 = 1;
	m_perm1[0] = 0.0; m_perm1[1] = 0.0; m_perm1[2] = 0.0;
	m_perm2[0] = 0.0; m_perm2[1] = 0.0; m_perm2[2] = 0.0;

	m_M0 = 0;
    m_alpha0 = 0;
	m_M[0] = 0.0; m_M[1] = 0.0; m_M[2] = 0.0;
	m_alpha[0] = 0.0; m_alpha[1] = 0.0; m_alpha[2] = 0.0;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
mat3ds FEPermRefOrtho::Permeability(FEMaterialPoint& mp)
{
	int a;
	vec3d V;			// orthonormal material directions in reference configuration
	mat3ds m[3];		// texture tensor in current configuration
	
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
    if (J <= phisr) feLogError("The perm-ref-ortho permeability calculation failed!\nThe volume ratio (J=%g) dropped below its theoretical minimum phi0=%g.",J,phisr);
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	for (a=0; a<3; a++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to V
		V.x = Q[0][a]; V.y = Q[1][a]; V.z = Q[2][a];
		m[a] = dyad(F*V);	// Evaluate texture tensor in the current configuration
	}
	
	// --- strain-dependent permeability ---
	
	double f;
    double M0 = m_M0(mp);
    double alpha0 = m_alpha0(mp);
	double k0 = m_perm0(mp)*pow((J-phisr)/(1-phi0),alpha0)*exp(M0*(J*J-1.0)/2.0);
	double k1[3] = { m_perm1[0](mp), m_perm1[1](mp), m_perm1[2](mp) };
	double k2[3] = { m_perm2[0](mp), m_perm2[1](mp), m_perm2[2](mp) };
	double alpha[3] = { m_alpha[0](mp), m_alpha[1](mp), m_alpha[2](mp) };
	double M[3] = { m_M[0](mp), m_M[1](mp), m_M[2](mp) };
	for (a=0; a<3; a++) {
		f = pow((J-phisr)/(1-phi0),alpha[a])*exp(M[a]*(J*J-1.0)/2.0);
		k1[a] *= f/(J*J);
		k2[a] *= 0.5*f/pow(J,4);
	}
	mat3ds kt = k0*I
	+k1[0]*m[0]+k1[1]*m[1]+k1[2]*m[2]
	+(2.0*k2[0])*(m[0]*b).sym()+(2.0*k2[1])*(m[1]*b).sym()+(2.0*k2[2])*(m[2]*b).sym();
	
	return kt;
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4dmm FEPermRefOrtho::Tangent_Permeability_Strain(FEMaterialPoint &mp)
{
	int a;
	vec3d V;			// orthonormal material directions in reference configuration
	mat3ds m[3];		// texture tensor in current configuration
	
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
    if (J <= phisr) feLogError("The perm-ref-ortho permeability calculation failed!\nThe volume ratio (J=%g) dropped below its theoretical minimum phi0=%g.",J,phisr);
    
	// get local coordinates
	mat3d Q = GetLocalCS(mp);

	for (a=0; a<3; a++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to V
		V.x = Q[0][a]; V.y = Q[1][a]; V.z = Q[2][a];
		m[a] = dyad(F*V);	// Evaluate texture tensor in the current configuration
	}
	
	double f, k0, K0prime, K1prime, K2prime;
	mat3ds k0hat, k1hat, k2hat;
    double M0 = m_M0(mp);
    double alpha0 = m_alpha0(mp);
	k0 = m_perm0(mp)*pow((J-phisr)/(1-phi0),alpha0)*exp(M0*(J*J-1.0)/2.0);
	K0prime = (1+J*(alpha0/(J-phisr)+M0*J))*k0;
	k0hat = mat3dd(K0prime);
	tens4dmm K4 = dyad1mm(I,k0hat)-dyad4s(I)*(2*k0);
	double k1[3] = { m_perm1[0](mp), m_perm1[1](mp), m_perm1[2](mp) };
	double k2[3] = { m_perm2[0](mp), m_perm2[1](mp), m_perm2[2](mp) };
	double alpha[3] = { m_alpha[0](mp), m_alpha[1](mp), m_alpha[2](mp) };
	double M[3] = { m_M[0](mp), m_M[1](mp), m_M[2](mp) };
	for (a=0; a<3; a++) {
		f = pow((J-phisr)/(1-phi0),alpha[a])*exp(M[a]*(J*J-1.0)/2.0);
		k1[a] *= f/(J*J);
		k2[a] *= 0.5*f/pow(J,4);
		K1prime = (J*J*M[a]+(J*(alpha[a]-1)+phi0)/(J-phisr))*k1[a];
		K2prime = (J*J*M[a]+(J*(alpha[a]-3)+3*phi0)/(J-phisr))*k2[a];
		k1hat = mat3dd(K1prime);
		k2hat = mat3dd(K2prime);
		K4 += dyad1mm(m[a],k1hat) + dyad1mm((m[a]*b).sym()*2.0,k2hat)
		+dyad4s(m[a],b)*(2.0*k2[a]);
	}
	
	return K4;
}
