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
#include "FECubicCLE.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FECubicCLE, FEElasticMaterial)
	ADD_PARAMETER(m_lp1, "lp1")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_lm1, "lm1")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_l2 , "l2")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_mu , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu")->setUnits(UNIT_PRESSURE);

    ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Check material parameters.
bool FECubicCLE::Validate()
{
    if (FEElasticMaterial::Validate() == false) return false;

	// Evaluate Lame coefficients
    double	lam[3][3];
    double	mu[3];
    mu[0] = m_mu;
    mu[1] = m_mu;
    mu[2] = m_mu;
    lam[0][0] = m_lm1; lam[0][1] = m_l2 ; lam[0][2] = m_l2;
    lam[1][0] = m_l2 ; lam[1][1] = m_lm1; lam[1][2] = m_l2;
    lam[2][0] = m_l2 ; lam[2][1] = m_l2 ; lam[2][2] = m_lm1;
    
    // check that stiffness matrix is positive definite
    mat3ds c(lam[0][0]+2*mu[0],lam[1][1]+2*mu[1],lam[2][2]+2*mu[2],lam[0][1],lam[1][2],lam[0][2]);
    double l[3];
    c.exact_eigen(l);
    
	if ((l[0] < 0) || (l[1] < 0) || (l[2] < 0))
	{
		feLogError("Stiffness matrix is not positive definite.");
		return false;
	}
    
    // repeat check with all tensile diagonal first lamé constants
    lam[0][0] = m_lp1; lam[1][1] = m_lp1; lam[2][2] = m_lp1;
    // check that compliance matrix is positive definite
    c = mat3ds(lam[0][0]+2*mu[0],lam[1][1]+2*mu[1],lam[2][2]+2*mu[2],lam[0][1],lam[1][2],lam[0][2]);
    c.exact_eigen(l);
    
	if ((l[0] < 0) || (l[1] < 0) || (l[2] < 0))
	{
		feLogError("Stiffness matrix is not positive definite.");
		return false;
	}
    
	return true;
}

//-----------------------------------------------------------------------------
//! Calculates the stress
mat3ds FECubicCLE::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // Evaluate Lame coefficients
    double	lam[3][3];
    double	mu[3];
    mu[0] = m_mu;
    mu[1] = m_mu;
    mu[2] = m_mu;
    lam[0][0] = m_lm1; lam[0][1] = m_l2 ; lam[0][2] = m_l2;
    lam[1][0] = m_l2 ; lam[1][1] = m_lm1; lam[1][2] = m_l2;
    lam[2][0] = m_l2 ; lam[2][1] = m_l2 ; lam[2][2] = m_lm1;
    
    int i,j;
    vec3d a0;           // texture direction in reference configuration
    mat3ds A0[3];		// texture tensor in reference configuration
    double AE[3];		// E:A
    mat3ds E = pt.Strain();
    
    mat3d F = pt.m_F;
    double J = pt.m_J;
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = Q[0][i]; a0.y = Q[1][i]; a0.z = Q[2][i];
        A0[i] = dyad(a0);
        AE[i] = a0*(E*a0);
    }
    
    lam[0][0] = (AE[0] >= 0) ? m_lp1 : m_lm1;
    lam[1][1] = (AE[1] >= 0) ? m_lp1 : m_lm1;
    lam[2][2] = (AE[2] >= 0) ? m_lp1 : m_lm1;
    
	mat3ds s; s.zero();
    
    for (i=0; i<3; ++i) {
        s += (A0[i]*E).sym()*(mu[i]);
        for (j=0; j<3; ++j) {
            s += A0[j]*(lam[i][j]*AE[i]);
        }
    }
    s = (F*s*F.transpose()).sym()/J;
    
    return s;
}

//-----------------------------------------------------------------------------
//! Calculates the elasticity tensor
tens4ds FECubicCLE::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // Evaluate Lame coefficients
    double	lam[3][3];
    double	mu[3];
    mu[0] = m_mu;
    mu[1] = m_mu;
    mu[2] = m_mu;
    lam[0][0] = m_lm1; lam[0][1] = m_l2 ; lam[0][2] = m_l2;
    lam[1][0] = m_l2 ; lam[1][1] = m_lm1; lam[1][2] = m_l2;
    lam[2][0] = m_l2 ; lam[2][1] = m_l2 ; lam[2][2] = m_lm1;
    
    int i,j;
    vec3d a0;           // texture direction in reference configuration
    mat3ds A[3];        // texture tensor in current configuration
    double AE[3];       // E:A
    mat3ds E = pt.Strain();

    mat3d F = pt.m_F;
    double J = pt.m_J;
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = Q[0][i]; a0.y = Q[1][i]; a0.z = Q[2][i];
        A[i] = dyad(F*a0);
        AE[i] = a0*(E*a0);
    }
    
    lam[0][0] = (AE[0] >= 0) ? m_lp1 : m_lm1;
    lam[1][1] = (AE[1] >= 0) ? m_lp1 : m_lm1;
    lam[2][2] = (AE[2] >= 0) ? m_lp1 : m_lm1;

    tens4ds c;
    c.zero();

    mat3ds B = pt.LeftCauchyGreen();
    for (i=0; i<3; ++i) {
        c += dyad4s(A[i], B)*mu[i];
        for (j=0; j<3; ++j) {
            c += dyad1s(A[i],A[j])*lam[i][j]/2;
        }
    }
    c /= J;
    
    return c;
}

//-----------------------------------------------------------------------------
//! Calculates the strain energy density
double FECubicCLE::StrainEnergyDensity(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // Evaluate Lame coefficients
    double	lam[3][3];
    double	mu[3];
    mu[0] = m_mu;
    mu[1] = m_mu;
    mu[2] = m_mu;
    lam[0][0] = m_lm1; lam[0][1] = m_l2 ; lam[0][2] = m_l2;
    lam[1][0] = m_l2 ; lam[1][1] = m_lm1; lam[1][2] = m_l2;
    lam[2][0] = m_l2 ; lam[2][1] = m_l2 ; lam[2][2] = m_lm1;
    
    int i,j;
    vec3d a0;           // texture direction in reference configuration
    mat3ds A0;          // texture tensor in reference configuration
    double AE[3], AE2[3];	// Ka
    mat3ds E = pt.Strain();
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = Q[0][i]; a0.y = Q[1][i]; a0.z = Q[2][i];
        A0 = dyad(a0);
        AE[i] = a0*(E*a0);
        AE2[i] = a0*(E*E*a0);
    }
    
    lam[0][0] = (AE[0] >= 0) ? m_lp1 : m_lm1;
    lam[1][1] = (AE[1] >= 0) ? m_lp1 : m_lm1;
    lam[2][2] = (AE[2] >= 0) ? m_lp1 : m_lm1;
    
    double W = 0;
    
    for (i=0; i<3; ++i) {
        W += AE2[i]*mu[i];
        for (j=0; j<3; ++j) {
            W += AE[i]*AE[j]*lam[i][j]/2;
        }
    }
    
    return W;
}
