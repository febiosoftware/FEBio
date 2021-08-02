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



#include "stdafx.h"
#include "FEPermRefOrtho.h"


// define the material parameters
BEGIN_FECORE_CLASS(FEPermRefOrtho, FEHydraulicPermeability)
	ADD_PARAMETER(m_perm0 , FE_RANGE_GREATER_OR_EQUAL(0.0), "perm0")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_M0    , FE_RANGE_GREATER_OR_EQUAL(0.0), "M0");
	ADD_PARAMETER(m_alpha0, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha0");
	ADD_PARAMETER(m_perm1, "perm1")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_perm2, "perm2")->setUnits(UNIT_PERMEABILITY);
	ADD_PARAMETER(m_M    , "M");
	ADD_PARAMETER(m_alpha, "alpha");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPermRefOrtho::FEPermRefOrtho(FEModel* pfem) : FEHydraulicPermeability(pfem)
{
	m_perm0 = 1;
	m_perm1 = vec3d(0,0,0);
    m_perm2 = vec3d(0,0,0);

	m_M0 = 0;
    m_alpha0 = 0;
	m_M = vec3d(0);
	m_alpha = vec3d(0);
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
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	for (a=0; a<3; a++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to V
		V.x = Q[0][a]; V.y = Q[1][a]; V.z = Q[2][a];
		m[a] = dyad(F*V);	// Evaluate texture tensor in the current configuration
	}
	
	// --- strain-dependent permeability ---
	
	double f, k1[3], k2[3], alpha[3], M[3];
    double M0 = m_M0(mp);
    double alpha0 = m_alpha0(mp);
	double k0 = m_perm0(mp)*pow((J-phi0)/(1-phi0),alpha0)*exp(M0*(J*J-1.0)/2.0);
    vec3d k1v = m_perm1(mp); k1[0] = k1v.x; k1[1] = k1v.y; k1[2] = k1v.z;
    vec3d k2v = m_perm2(mp); k2[0] = k2v.x; k2[1] = k2v.y; k2[2] = k2v.z;
    vec3d alphav = m_alpha(mp); alpha[0] = alphav.x; alpha[1] = alphav.y; alpha[2] = alphav.z;
    vec3d Mv = m_M(mp); M[0] = Mv.x; M[1] = Mv.y; M[2] = Mv.z;
	for (a=0; a<3; a++) {
		f = pow((J-phi0)/(1-phi0),alpha[a])*exp(M[a]*(J*J-1.0)/2.0);
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
	
	// get local coordinates
	mat3d Q = GetLocalCS(mp);

	for (a=0; a<3; a++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to V
		V.x = Q[0][a]; V.y = Q[1][a]; V.z = Q[2][a];
		m[a] = dyad(F*V);	// Evaluate texture tensor in the current configuration
	}
	
	double f, k0, k1[3], k2[3], alpha[3], M[3], K0prime, K1prime, K2prime;
	mat3ds k0hat, k1hat, k2hat;
    double M0 = m_M0(mp);
    double alpha0 = m_alpha0(mp);
	k0 = m_perm0(mp)*pow((J-phi0)/(1-phi0),alpha0)*exp(M0*(J*J-1.0)/2.0);
	K0prime = (1+J*(alpha0/(J-phi0)+M0*J))*k0;
	k0hat = mat3dd(K0prime);
	tens4dmm K4 = dyad1mm(I,k0hat)-dyad4s(I)*(2*k0);
    vec3d k1v = m_perm1(mp); k1[0] = k1v.x; k1[1] = k1v.y; k1[2] = k1v.z;
    vec3d k2v = m_perm2(mp); k2[0] = k2v.x; k2[1] = k2v.y; k2[2] = k2v.z;
    vec3d alphav = m_alpha(mp); alpha[0] = alphav.x; alpha[1] = alphav.y; alpha[2] = alphav.z;
    vec3d Mv = m_M(mp); M[0] = Mv.x; M[1] = Mv.y; M[2] = Mv.z;
	for (a=0; a<3; a++) {
		f = pow((J-phi0)/(1-phi0),alpha[a])*exp(M[a]*(J*J-1.0)/2.0);
		k1[a] *= f/(J*J);
		k2[a] *= 0.5*f/pow(J,4);
		K1prime = (J*J*M[a]+(J*(alpha[a]-1)+phi0)/(J-phi0))*k1[a];
		K2prime = (J*J*M[a]+(J*(alpha[a]-3)+3*phi0)/(J-phi0))*k2[a];
		k1hat = mat3dd(K1prime);
		k2hat = mat3dd(K2prime);
		K4 += dyad1mm(m[a],k1hat) + dyad1mm((m[a]*b).sym()*2.0,k2hat)
		+dyad4s(m[a],b)*(2.0*k2[a]);
	}
	
	return K4;
}
