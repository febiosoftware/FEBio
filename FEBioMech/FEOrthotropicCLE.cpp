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
#include "FEOrthotropicCLE.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEOrthotropicCLE, FEElasticMaterial)
	ADD_PARAMETER(lp11, "lp11");
	ADD_PARAMETER(lm11, "lm11");
	ADD_PARAMETER(lp22, "lp22");
	ADD_PARAMETER(lm22, "lm22");
	ADD_PARAMETER(lp33, "lp33");
	ADD_PARAMETER(lm33, "lm33");
	ADD_PARAMETER(l12, "l12");
	ADD_PARAMETER(l23, "l23");
	ADD_PARAMETER(l31, "l31");
	ADD_PARAMETER(mu1, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu1");
	ADD_PARAMETER(mu2, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu2");
	ADD_PARAMETER(mu3, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu3");

    ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Check material parameters.
bool FEOrthotropicCLE::Validate()
{
    if (FEElasticMaterial::Validate() == false) return false;
    
    // Evaluate Lame coefficients
    double	lam[3][3];
    double	mu[3];
    mu[0] = mu1;
    mu[1] = mu2;
    mu[2] = mu3;
    lam[0][0] = lm11; lam[0][1] = l12 ; lam[0][2] = l31;
    lam[1][0] = l12 ; lam[1][1] = lm22; lam[1][2] = l23;
    lam[2][0] = l31 ; lam[2][1] = l23 ; lam[2][2] = lm33;
    
    // check that stiffness matrix is positive definite
    mat3ds c(lam[0][0]+2*mu[0],lam[1][1]+2*mu[1],lam[2][2]+2*mu[2],lam[0][1],lam[1][2],lam[0][2]);
    double l[3];
    c.exact_eigen(l);
    
	if ((l[0] < 0) || (l[1] < 0) || (l[2] < 0)) {
		feLogError("Stiffness matrix is not positive definite.");
		return false;
	}
    
    // repeat check with all tensile diagonal first lamé constants
    lam[0][0] = lp11; lam[1][1] = lp22; lam[2][2] = lp33;
    // check that compliance matrix is positive definite
    c = mat3ds(lam[0][0]+2*mu[0],lam[1][1]+2*mu[1],lam[2][2]+2*mu[2],lam[0][1],lam[1][2],lam[0][2]);
    c.exact_eigen(l);
    
	if ((l[0] < 0) || (l[1] < 0) || (l[2] < 0)) {
		feLogError("Stiffness matrix is not positive definite.");
		return false;
	}
 
	return true;
}

//-----------------------------------------------------------------------------
//! Calculates the stress
mat3ds FEOrthotropicCLE::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // Evaluate Lame coefficients
    double	lam[3][3];
    double	mu[3];
    mu[0] = mu1;
    mu[1] = mu2;
    mu[2] = mu3;
    lam[0][0] = lm11; lam[0][1] = l12 ; lam[0][2] = l31;
    lam[1][0] = l12 ; lam[1][1] = lm22; lam[1][2] = l23;
    lam[2][0] = l31 ; lam[2][1] = l23 ; lam[2][2] = lm33;
    
    int i,j;
    vec3d a0;           // texture direction in reference configuration
    mat3ds A[3];		// texture tensor in current configuration
    double K[3];		// Ka
    mat3ds BmI = pt.LeftCauchyGreen() - mat3dd(1);			// B - I
    
    mat3d F = pt.m_F;
    double J = pt.m_J;

	// get local coordinates
	mat3d Q = GetLocalCS(mp);

    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = Q[0][i]; a0.y = Q[1][i]; a0.z = Q[2][i];
        A[i] = dyad(F*a0);
        K[i] = 0.5*(A[i].tr() - 1);
    }
    
    lam[0][0] = (K[0] >= 0) ? lp11 : lm11;
    lam[1][1] = (K[1] >= 0) ? lp22 : lm22;
    lam[2][2] = (K[2] >= 0) ? lp33 : lm33;
    
	mat3ds s; 
	s.zero();
    
    for (i=0; i<3; ++i) {
        s += (A[i]*BmI).sym()*(mu[i]);
        for (j=0; j<3; ++j) {
            s += A[j]*(lam[i][j]*K[i]);
        }
    }
    s /= J;
    
    return s;
}

//-----------------------------------------------------------------------------
//! Calculates the elasticity tensor
tens4ds FEOrthotropicCLE::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // Evaluate Lame coefficients
    double	lam[3][3];
    double	mu[3];
    mu[0] = mu1;
    mu[1] = mu2;
    mu[2] = mu3;
    lam[0][0] = lm11; lam[0][1] = l12 ; lam[0][2] = l31;
    lam[1][0] = l12 ; lam[1][1] = lm22; lam[1][2] = l23;
    lam[2][0] = l31 ; lam[2][1] = l23 ; lam[2][2] = lm33;
    
    int i,j;
    vec3d a0;           // texture direction in reference configuration
    mat3ds A[3];		// texture tensor in current configuration
    double K[3];		// Ka
    mat3ds B = pt.LeftCauchyGreen();
    
    mat3d F = pt.m_F;
    double J = pt.m_J;
    
	// get local coordinates
	mat3d Q = GetLocalCS(mp);

    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = Q[0][i]; a0.y = Q[1][i]; a0.z = Q[2][i];
        A[i] = dyad(F*a0);
        K[i] = 0.5*(A[i].tr() - 1);
    }
    
    lam[0][0] = (K[0] >= 0) ? lp11 : lm11;
    lam[1][1] = (K[1] >= 0) ? lp22 : lm22;
    lam[2][2] = (K[2] >= 0) ? lp33 : lm33;
    
    tens4ds c;
    c.zero();
    
    for (i=0; i<3; ++i) {
        c += dyad4s(A[i], B)*mu[i];
        for (j=0; j<3; ++j) {
            c += dyad1s(A[i],A[j])*lam[i][j];
        }
    }
    c /= J;
    
    return c;
}

//-----------------------------------------------------------------------------
//! Calculates the strain energy density
double FEOrthotropicCLE::StrainEnergyDensity(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // Evaluate Lame coefficients
    double	lam[3][3];
    double	mu[3];
    mu[0] = mu1;
    mu[1] = mu2;
    mu[2] = mu3;
    lam[0][0] = lm11; lam[0][1] = l12 ; lam[0][2] = l31;
    lam[1][0] = l12 ; lam[1][1] = lm22; lam[1][2] = l23;
    lam[2][0] = l31 ; lam[2][1] = l23 ; lam[2][2] = lm33;
    
    int i,j;
    vec3d a0;           // texture direction in reference configuration
    mat3ds A0;          // texture tensor in current configuration
    double K[3], L[3];	// Ka
    mat3ds E = (pt.RightCauchyGreen() - mat3dd(1))/2;
    
	// get local coordinates
	mat3d Q = GetLocalCS(mp);

    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = Q[0][i]; a0.y = Q[1][i]; a0.z = Q[2][i];
        A0 = dyad(a0);
        K[i] = (A0*E).trace();
        L[i] = (A0*E*E).trace();
    }
    
    lam[0][0] = (K[0] >= 0) ? lp11 : lm11;
    lam[1][1] = (K[1] >= 0) ? lp22 : lm22;
    lam[2][2] = (K[2] >= 0) ? lp33 : lm33;
    
    double W = 0;
    
    for (i=0; i<3; ++i) {
        W += L[i]*mu[i];
        for (j=0; j<3; ++j) {
            W += K[i]*K[j]*lam[i][j]/2;
        }
    }
    
    return W;
}
