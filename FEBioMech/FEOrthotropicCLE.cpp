//
//  FEOrthotropicCLE.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/25/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FEOrthotropicCLE.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEOrthotropicCLE, FEElasticMaterial)
	ADD_PARAMETER(lp11, FE_PARAM_DOUBLE, "lp11");
	ADD_PARAMETER(lm11, FE_PARAM_DOUBLE, "lm11");
	ADD_PARAMETER(lp22, FE_PARAM_DOUBLE, "lp22");
	ADD_PARAMETER(lm22, FE_PARAM_DOUBLE, "lm22");
	ADD_PARAMETER(lp33, FE_PARAM_DOUBLE, "lp33");
	ADD_PARAMETER(lp11, FE_PARAM_DOUBLE, "lp11");
	ADD_PARAMETER(lm33, FE_PARAM_DOUBLE, "lm33");
	ADD_PARAMETER(l12, FE_PARAM_DOUBLE, "l12");
	ADD_PARAMETER(l23, FE_PARAM_DOUBLE, "l23");
	ADD_PARAMETER(l31, FE_PARAM_DOUBLE, "l31");
	ADD_PARAMETER2(mu1, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu1");
	ADD_PARAMETER2(mu2, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu2");
	ADD_PARAMETER2(mu3, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu3");
END_PARAMETER_LIST();

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
    
    if ((l[0]<0) || (l[1]<0) || (l[2]<0))
        return MaterialError("Stiffness matrix is not positive definite.");
    
    // repeat check with all tensile diagonal first lamé constants
    lam[0][0] = lp11; lam[1][1] = lp22; lam[2][2] = lp33;
    // check that compliance matrix is positive definite
    c = mat3ds(lam[0][0]+2*mu[0],lam[1][1]+2*mu[1],lam[2][2]+2*mu[2],lam[0][1],lam[1][2],lam[0][2]);
    c.exact_eigen(l);
    
    if ((l[0]<0) || (l[1]<0) || (l[2]<0))
        return MaterialError("Stiffness matrix is not positive definite.");
 
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
    
    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = pt.m_Q[0][i]; a0.y = pt.m_Q[1][i]; a0.z = pt.m_Q[2][i];
        A[i] = dyad(F*a0);
        K[i] = 0.5*(A[i].tr() - 1);
    }
    
    lam[0][0] = (K[0] >= 0) ? lp11 : lm11;
    lam[1][1] = (K[1] >= 0) ? lp22 : lm22;
    lam[2][2] = (K[2] >= 0) ? lp33 : lm33;
    
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
    
    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = pt.m_Q[0][i]; a0.y = pt.m_Q[1][i]; a0.z = pt.m_Q[2][i];
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
    
    for (i=0; i<3; i++) {	// Perform sum over all three texture directions
        // Copy the texture direction in the reference configuration to a0
        a0.x = pt.m_Q[0][i]; a0.y = pt.m_Q[1][i]; a0.z = pt.m_Q[2][i];
        A0 = dyad(a0);
        K[i] = (A0*E).tr();
        L[i] = (A0*E*E).tr();
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
