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
#include "FEOrthoElastic.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEOrthoElastic, FEElasticMaterial)
	ADD_PARAMETER(m_E1 , FE_RANGE_GREATER(0.0), "E1")->setLongName("E1 modulus");
	ADD_PARAMETER(m_E2 , FE_RANGE_GREATER(0.0), "E2")->setLongName("E2 modulus");
	ADD_PARAMETER(m_E3 , FE_RANGE_GREATER(0.0), "E3")->setLongName("E3 modulus");
	ADD_PARAMETER(m_G12, FE_RANGE_GREATER_OR_EQUAL(0.0), "G12")->setLongName("G12 shear modulus");
	ADD_PARAMETER(m_G23, FE_RANGE_GREATER_OR_EQUAL(0.0), "G23")->setLongName("G23 shear modulus");
	ADD_PARAMETER(m_G31, FE_RANGE_GREATER_OR_EQUAL(0.0), "G31")->setLongName("G31 shear modulus");
	ADD_PARAMETER(m_v12, "v12");
	ADD_PARAMETER(m_v23, "v23");
	ADD_PARAMETER(m_v31, "v31");

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEOrthoElastic::FEOrthoElastic(FEModel* pfem) : FEElasticMaterial(pfem) {}

//-----------------------------------------------------------------------------
//! Material initialization.
bool FEOrthoElastic::Validate()
{
	if (FEElasticMaterial::Validate() == false) return false;

	// do some sanity checks
	if (m_v12.isConst() && m_E1.isConst() && m_E2.isConst())
	{
		double v12 = m_v12.constValue();
		double E1 = m_E1.constValue();
		double E2 = m_E2.constValue();
		if (v12 > sqrt(E1 / E2)) { feLogError("Invalid value for v12. Let v12 <= sqrt(E1/E2)"); return false; }
	}

	if (m_v23.isConst() && m_E2.isConst() && m_E3.isConst())
	{
		double v23 = m_v23.constValue();
		double E2 = m_E2.constValue();
		double E3 = m_E3.constValue();

		if (v23 > sqrt(E2 / E3)) { feLogError("Invalid value for v23. Let v23 <= sqrt(E2/E3)"); return false; }
	}

	if (m_v31.isConst() && m_E1.isConst() && m_E3.isConst())
	{
		double v31 = m_v31.constValue();
		double E1 = m_E1.constValue();
		double E3 = m_E3.constValue();

		if (v31 > sqrt(E3 / E1)) { feLogError("Invalid value for v31. Let v31 <= sqrt(E3/E1)"); return false; }
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEOrthoElastic::EvaluateLameCoefficients(FEMaterialPoint& mp, double lam[3][3], double mu[3])
{
	double E1  = m_E1 (mp); double  E2 = m_E2 (mp); double E3  = m_E3 (mp);
	double v12 = m_v12(mp); double v23 = m_v23(mp); double v31 = m_v31(mp);
	double G12 = m_G12(mp); double G23 = m_G23(mp); double G31 = m_G31(mp);

	// Evaluate Lame coefficients
	mu[0] = G12 + G31 - G23;
	mu[1] = G12 - G31 + G23;
	mu[2] =-G12 + G31 + G23;
	lam[0][0] = 1.0/E1; lam[0][1] = -v12/E1; lam[0][2] = -v31/E3;
	lam[1][0] = -v12/E1; lam[1][1] = 1.0/E2; lam[1][2] = -v23/E2;
	lam[2][0] = -v31/E3; lam[2][1] = -v23/E2; lam[2][2] = 1.0/E3;

	// check that compliance matrix is positive definite
	mat3ds c(lam[0][0], lam[1][1], lam[2][2], lam[0][1], lam[1][2], lam[0][2]);
	double l[3];
	c.exact_eigen(l);

	if ((l[0] < 0) || (l[1] < 0) || (l[2] < 0)) {
		feLogError("Stiffness matrix is not positive definite.");
		assert(false);
		return false;
	}

	// evaluate stiffness matrix and extract Lame constants
	c = c.inverse();
	lam[0][0] = c(0,0) - 2*mu[0];
	lam[1][1] = c(1,1) - 2*mu[1];
	lam[2][2] = c(2,2) - 2*mu[2];
	lam[1][2] = c(1,2); lam[2][1] = c(2,1);
	lam[2][0] = c(2,0); lam[0][2] = c(0,2);
	lam[0][1] = c(0,1); lam[1][0] = c(1,0);

	return true;
}

//-----------------------------------------------------------------------------
//! Calculates the stress for a linear orthotropic material
mat3ds FEOrthoElastic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	int i,j;
	vec3d a0[3];		// texture direction in reference configuration
	vec3d a[3];			// texture direction in current configuration
	mat3ds A[3];		// texture tensor in current configuration
	double K[3];		// Ka
	double L[3];		// La
	mat3ds bmi;			// B - I
	// Evaluate the deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// calculate left and right Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	mat3ds c = pt.RightCauchyGreen();
	mat3ds c2 = c.sqr();
	mat3dd I(1.);
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	for (i=0; i<3; i++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to a0
		a0[i].x = Q[0][i]; a0[i].y = Q[1][i]; a0[i].z = Q[2][i];
		K[i] = a0[i]*(c*a0[i]);
		L[i] = a0[i]*(c2*a0[i]);
		a[i] = F*a0[i]/sqrt(K[i]);	// Evaluate the texture direction in the current configuration
		A[i] = dyad(a[i]);			// Evaluate the texture tensor in the current configuration
	}

	// evaluate the Lame coefficients
	double lam[3][3] = { 0 };
	double mu[3] = { 0 };
	EvaluateLameCoefficients(mp, lam, mu);
	
	// Evaluate the stress
	mat3ds s;
	s.zero();		// Initialize for summation
	bmi = b - I;
	for (i=0; i<3; i++) {
		s += (A[i]*bmi).sym()*(2.0*mu[i] * K[i]);
		for (j=0; j<3; j++)
			s += lam[i][j]*((K[i]-1)*K[j]*A[j]+(K[j]-1)*K[i]*A[i])/2.;
	}
	s /= 2.0*J;
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEOrthoElastic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	int i,j;
	vec3d a0[3];		// texture direction in reference configuration
	vec3d a[3];			// texture direction in current configuration
	mat3ds A[3];		// texture tensor in current configuration
	double K[3];		// Ka
	// Evaluate the strain and texture
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// calculate left and right Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	mat3ds c = pt.RightCauchyGreen();
	mat3dd I(1.);
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	for (i=0; i<3; i++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to a0
		a0[i].x = Q[0][i]; a0[i].y = Q[1][i]; a0[i].z = Q[2][i];
		K[i] = a0[i]*(c*a0[i]);
		a[i] = F*a0[i]/sqrt(K[i]);	// Evaluate the texture direction in the current configuration
		A[i] = dyad(a[i]);			// Evaluate the texture tensor in the current configuration
	}
	
	// evaluate the Lame coefficients
	double lam[3][3] = { 0 };
	double mu[3] = { 0 };
	EvaluateLameCoefficients(mp, lam, mu);

	tens4ds C(0.0);
	for (i=0; i<3; i++) {
		C += mu[i]*K[i]*dyad4s(A[i],b);
		for (j=0; j<3; j++)
			C += lam[i][j]*K[i]*K[j]*dyad1s(A[i],A[j])/2.;
	}
	
	// Elasticity tensor
	C /= J;
	
	return C;
}

//-----------------------------------------------------------------------------
double FEOrthoElastic::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    mat3ds E = (pt.RightCauchyGreen() - mat3dd(1))/2;
    mat3ds E2 = E.sqr();
    
	vec3d a0[3];		// texture direction in reference configuration
	mat3ds A0[3];		// texture tensor in current configuration
    double AE[3], AE2[3];

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	for (int i=0; i<3; i++) {
		// Copy the texture direction in the reference configuration to a0
		a0[i].x = Q[0][i]; a0[i].y = Q[1][i]; a0[i].z = Q[2][i];
		A0[i] = dyad(a0[i]);			// Evaluate the texture tensor in the reference configuration
        AE[i] = A0[i].dotdot(E);
        AE2[i] = A0[i].dotdot(E2);
	}
    
	// evaluate the Lame coefficients
	double lam[3][3] = { 0 };
	double mu[3] = { 0 };
	EvaluateLameCoefficients(mp, lam, mu);

	// calculate strain energy
    double sed = mu[0]*AE2[0] + mu[1]*AE2[1] + mu[2]*AE2[2]
    +0.5*(lam[0][0]*AE[0]*AE[0]+lam[1][1]*AE[1]*AE[1]+lam[2][2]*AE[2]*AE[2])
    +lam[0][1]*AE[0]*AE[1]+lam[1][2]*AE[1]*AE[2]+lam[2][0]*AE[2]*AE[0];
    
	return sed;
}
