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
#include "FEFungOrthotropic.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEFungOrthotropic, FEUncoupledMaterial)
	ADD_PARAMETER(E1, FE_RANGE_GREATER(0.0), "E1")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(E2, FE_RANGE_GREATER(0.0), "E2")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(E3, FE_RANGE_GREATER(0.0), "E3")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(G12, FE_RANGE_GREATER_OR_EQUAL(0.0), "G12")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(G23, FE_RANGE_GREATER_OR_EQUAL(0.0), "G23")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(G31, FE_RANGE_GREATER_OR_EQUAL(0.0), "G31")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(v12, "v12");
	ADD_PARAMETER(v23, "v23");
	ADD_PARAMETER(v31, "v31");
	ADD_PARAMETER(m_c, FE_RANGE_GREATER(0.0), "c")->setUnits(UNIT_PRESSURE);

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFungOrthotropic::FEFungOrthotropic(FEModel* pfem) : FEUncoupledMaterial(pfem) 
{
}

//-----------------------------------------------------------------------------
//! Data initialization
bool FEFungOrthotropic::Validate()
{
	if (FEUncoupledMaterial::Validate() == false) return false;

	if (v12 > sqrt(E1/E2)) { feLogError("Invalid value for v12. Let v12 <= sqrt(E1/E2)"); return false; }
	if (v23 > sqrt(E2/E3)) { feLogError("Invalid value for v23. Let v23 <= sqrt(E2/E3)"); return false; }
	if (v31 > sqrt(E3/E1)) { feLogError("Invalid value for v31. Let v31 <= sqrt(E3/E1)"); return false; }
	
	// Evaluate Lame coefficients
	mu[0] = G12 + G31 - G23;
	mu[1] = G12 - G31 + G23;
	mu[2] =-G12 + G31 + G23;
	lam[0][0] = 1.0/E1; lam[0][1] = -v12/E1; lam[0][2] = -v31/E3;
	lam[1][0] = -v12/E1; lam[1][1] = 1.0/E2; lam[1][2] = -v23/E2;
	lam[2][0] = -v31/E3; lam[2][1] = -v23/E2; lam[2][2] = 1.0/E3;
	mat3d c(lam);
	c = c.inverse();
	lam[0][0] = c[0][0] - 2*mu[0];
	lam[1][1] = c[1][1] - 2*mu[1];
	lam[2][2] = c[2][2] - 2*mu[2];
	lam[1][2] = c[1][2]; lam[2][1] = c[2][1];
	lam[2][0] = c[2][0]; lam[0][2] = c[0][2];
	lam[0][1] = c[0][1]; lam[1][0] = c[1][0];

	return true;
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric stress
mat3ds FEFungOrthotropic::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	int i,j;
	vec3d a0[3];		// texture direction in reference configuration
	vec3d a[3];			// texture direction in current configuration
	mat3ds A[3];		// texture tensor in current configuration
	double K[3];		// Ka
	double L[3];		// La
	double eQ;			// exp(Q)
	mat3ds bmi;			// B - I
	// Evaluate the deformation gradient
	mat3d &F = pt.m_F;
	double detF = pt.m_J;
	mat3d Fd = F*pow(detF,-1./3.);
	
	// calculate deviatoric left and right Cauchy-Green tensor
	mat3ds b = pt.DevLeftCauchyGreen();
	mat3ds c = pt.DevRightCauchyGreen();
	mat3ds c2 = c.sqr();
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	for (i=0; i<3; i++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to a0
		a0[i].x = Q[0][i]; a0[i].y = Q[1][i]; a0[i].z = Q[2][i];
		K[i] = a0[i]*(c*a0[i]);
		L[i] = a0[i]*(c2*a0[i]);
		a[i] = Fd*a0[i]/sqrt(K[i]);	// Evaluate the texture direction in the current configuration
		A[i] = dyad(a[i]);			// Evaluate the texture tensor in the current configuration
	}
	
	// Evaluate exp(Q)
	eQ = 0;
	for (i=0; i<3; i++) {
		eQ += 2*mu[i]*(L[i]-2*K[i]+1);
		for (j=0; j<3; j++)
			eQ += lam[i][j]*(K[i]-1)*(K[j]-1);
	}
	eQ = exp(eQ/(4.*m_c));
	
	// Evaluate the stress
	mat3ds s;
	s.zero();		// Initialize for summation
	bmi = b - mat3dd(1.);
	for (i=0; i<3; i++) {
		s += (A[i]*bmi).sym()*(2.0*mu[i] * K[i]);
		for (j=0; j<3; j++)
			s += lam[i][j]*((K[i]-1)*K[j]*A[j]+(K[j]-1)*K[i]*A[i])/2.;
	}
	s *= eQ/(2.0*detF);

	return s.dev();
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric tangent
tens4ds FEFungOrthotropic::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	int i,j;
	vec3d a0[3];		// texture direction in reference configuration
	vec3d a[3];			// texture direction in current configuration
	mat3ds A[3];		// texture tensor in current configuration
	double K[3];		// Ka
	double L[3];		// La
	double eQ;			// exp(Q)
	mat3ds bmi;			// B - I
	// Evaluate the strain and texture
	mat3d &F = pt.m_F;
	double detF = pt.m_J;
	mat3d Fd = F*pow(detF,-1./3.);
	
	// calculate left and right Cauchy-Green tensor
	mat3ds b = pt.DevLeftCauchyGreen();
	mat3ds c = pt.DevRightCauchyGreen();
	mat3ds c2 = c.sqr();
	mat3dd I(1.);
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	for (i=0; i<3; i++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to a0
		a0[i].x = Q[0][i]; a0[i].y = Q[1][i]; a0[i].z = Q[2][i];
		K[i] = a0[i]*(c*a0[i]);
		L[i] = a0[i]*(c2*a0[i]);
		a[i] = Fd*a0[i]/sqrt(K[i]);	// Evaluate the texture direction in the current configuration
		A[i] = dyad(a[i]);			// Evaluate the texture tensor in the current configuration
	}
	
	// Evaluate exp(Q)
	eQ = 0;
	for (i=0; i<3; i++) {
		eQ += 2*mu[i]*(L[i]-2*K[i]+1);
		for (j=0; j<3; j++)
			eQ += lam[i][j]*(K[i]-1)*(K[j]-1);
	}
	eQ = exp(eQ/(4.*m_c));
	
	// Evaluate the distortional part of the Cauchy stress
	mat3ds sd;
	sd.zero();		// Initialize for summation
	tens4ds C(0.0);
	bmi = b - I;
	for (i=0; i<3; i++) {
		sd += ((A[i]*bmi).sym())*(2.0*mu[i] * K[i]);
		for (j=0; j<3; j++)
			sd += lam[i][j]*((K[i]-1)*K[j]*A[j]+(K[j]-1)*K[i]*A[i])/2.;
		C += mu[i]*K[i]*dyad4s(A[i],b);
		for (j=0; j<3; j++)
			C += lam[i][j]*K[i]*K[j]*dyad1s(A[i],A[j])/2.;
	}
	sd *= eQ/(2.0*detF);
	// This is the distortional part of the elasticity tensor
	C = (eQ/detF)*C + (2*detF/m_c/eQ)*dyad1s(sd);
	// This is the final value of the elasticity tensor
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);

	C += - 1./3.*(ddots(C,IxI) - IxI*(C.tr()/3.))
	+ 2./3.*((I4-IxI/3.)*sd.tr()-dyad1s(sd.dev(),I));
	
	return C;
}

//-----------------------------------------------------------------------------
//! Calculates the strain energy density
double FEFungOrthotropic::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	int i,j;
	vec3d a0[3];		// texture direction in reference configuration
	mat3ds A0[3];		// texture tensor in current configuration
	double eQ = 0.0;	// exp(Q)
    double AE[3], AE2[3];
	
	// calculate right Cauchy-Green tensor and Lagrange strain tensor
	mat3ds C = pt.DevRightCauchyGreen();
	mat3dd I(1.);
    mat3ds E = (C - I)*0.5;
    mat3ds E2 = E.sqr();
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	for (i=0; i<3; i++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to a0
		a0[i].x = Q[0][i]; a0[i].y = Q[1][i]; a0[i].z = Q[2][i];
		A0[i] = dyad(a0[i]);			// Evaluate the texture tensor in the reference configuration
        AE[i] = A0[i].dotdot(E);
        AE2[i] = A0[i].dotdot(E2);
	}
	
	// Evaluate exp(Q)
	for (i=0; i<3; i++) {
		eQ += 2*mu[i]*AE2[i];
		for (j=0; j<3; j++)
			eQ += lam[i][j]*AE[i]*AE[j];
	}
	eQ = exp(eQ/m_c);
	
	// Evaluate the strain energy density
	double sed = 0.5*m_c*(eQ-1);
	
	return sed;
}
