//
//  FECubicCLE.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/26/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FECubicCLE.h"
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FECubicCLE, FEElasticMaterial)
	ADD_PARAMETER(m_lp1, FE_PARAM_DOUBLE, "lp1");
	ADD_PARAMETER(m_lm1, FE_PARAM_DOUBLE, "lm1");
	ADD_PARAMETER(m_l2 , FE_PARAM_DOUBLE, "l2");
	ADD_PARAMETER2(m_mu , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Check material parameters.
void FECubicCLE::Init()
{
    FEElasticMaterial::Init();

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
    
    if ((l[0]<0) || (l[1]<0) || (l[2]<0))
        throw MaterialError("Stiffness matrix is not positive definite.");
    
    // repeat check with all tensile diagonal first lamé constants
    lam[0][0] = m_lp1; lam[1][1] = m_lp1; lam[2][2] = m_lp1;
    // check that compliance matrix is positive definite
    c = mat3ds(lam[0][0]+2*mu[0],lam[1][1]+2*mu[1],lam[2][2]+2*mu[2],lam[0][1],lam[1][2],lam[0][2]);
    c.exact_eigen(l);
    
    if ((l[0]<0) || (l[1]<0) || (l[2]<0))
        throw MaterialError("Stiffness matrix is not positive definite.");
    
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
    mat3ds A[3];		// texture tensor in current configuration
    double K[3];		// Ka
    mat3ds BmI = pt.LeftCauchyGreen() - mat3dd(1);			// B - I
    
    mat3d F = pt.m_F;
    double J = pt.m_J;
    
    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = pt.m_Q[0][i]; a0.y = pt.m_Q[1][i]; a0.z = pt.m_Q[2][i];
        A[i] = dyad(F*a0);
        K[i] = 0.5*(A[i].tr() - 1);
    }
    
    lam[0][0] = (K[0] >= 0) ? m_lp1 : m_lm1;
    lam[1][1] = (K[1] >= 0) ? m_lp1 : m_lm1;
    lam[2][2] = (K[2] >= 0) ? m_lp1 : m_lm1;
    
    mat3ds s(0);
    
    for (i=0; i<3; ++i) {
        s += (A[i]*BmI+BmI*A[i])*(mu[i]/2);
        for (j=0; j<3; ++j) {
            s += A[j]*(lam[i][j]*K[i]);
        }
    }
    s /= J;
    
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
    mat3ds A[3];		// texture tensor in current configuration
    double K[3];		// Ka
    mat3ds B = pt.LeftCauchyGreen();
    
    mat3d F = pt.m_F;
    double J = pt.m_J;
    
    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = pt.m_Q[0][i]; a0.y = pt.m_Q[1][i]; a0.z = pt.m_Q[2][i];
        A[i] = dyad(F*a0);
        K[i] = 0.5*(A[i].tr() - 1);
    }
    
    lam[0][0] = (K[0] >= 0) ? m_lp1 : m_lm1;
    lam[1][1] = (K[1] >= 0) ? m_lp1 : m_lm1;
    lam[2][2] = (K[2] >= 0) ? m_lp1 : m_lm1;
    
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
    mat3ds A0;          // texture tensor in current configuration
    double K[3], L[3];	// Ka
    mat3ds E = (pt.RightCauchyGreen() - mat3dd(1))/2;
    
    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = pt.m_Q[0][i]; a0.y = pt.m_Q[1][i]; a0.z = pt.m_Q[2][i];
        A0 = dyad(a0);
        K[i] = (A0*E).tr();
        L[i] = (A0*E*E).tr();
    }
    
    lam[0][0] = (K[0] >= 0) ? m_lp1 : m_lm1;
    lam[1][1] = (K[1] >= 0) ? m_lp1 : m_lm1;
    lam[2][2] = (K[2] >= 0) ? m_lp1 : m_lm1;
    
    double W = 0;
    
    for (i=0; i<3; ++i) {
        W += L[i]*mu[i];
        for (j=0; j<3; ++j) {
            W += K[i]*K[j]*lam[i][j]/2;
        }
    }
    
    return W;
}
