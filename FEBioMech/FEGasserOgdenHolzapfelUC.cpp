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
#include "FEGasserOgdenHolzapfelUC.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEGasserOgdenHolzapfelUC, FEUncoupledMaterial)
	ADD_PARAMETER(m_c    , FE_RANGE_GREATER_OR_EQUAL(0.0), "c");
	ADD_PARAMETER(m_k1   , FE_RANGE_GREATER_OR_EQUAL(0.0), "k1");
	ADD_PARAMETER(m_k2   , FE_RANGE_GREATER_OR_EQUAL(0.0), "k2");
	ADD_PARAMETER(m_kappa, FE_RANGE_CLOSED(0.0, 1.0/3.0), "kappa");
	ADD_PARAMETER(m_g    , "gamma");
END_FECORE_CLASS();

#define ONE 0.9999

//-----------------------------------------------------------------------------
//! Calculates the deviatoric stress
mat3ds FEGasserOgdenHolzapfelUC::DevStress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // determinant of deformation gradient
    double J = pt.m_J;

   
    // Evaluate the distortional deformation gradient
	double Jm13 = pow(J, -1. / 3.);
	mat3d F = pt.m_F*Jm13;
    
    // calculate deviatoric left Cauchy-Green tensor: b = F*Ft
    mat3ds b = pt.LeftCauchyGreen()*(Jm13*Jm13);
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    // Copy the local element basis directions to n
	vec3d n[2];
    n[0].x = Q[0][0]; n[0].y = Q[1][0]; n[0].z = Q[2][0];
    n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1];
    
    // Evaluate the structural direction in the current configuration
    double cg = cos(m_g); double sg = sin(m_g);
	vec3d a[2];
    a[0] = F*(n[0]*cg + n[1]*sg);
    a[1] = F*(n[0]*cg - n[1]*sg);
    
    // Evaluate the ground matrix stress
    mat3ds s = m_c*b;
    
    // Evaluate the structural tensors in the current configuration
    // and the fiber strains and stress contributions
    mat3ds h0 = m_kappa*b;
    if (a[0]*a[0] > ONE)
        h0 += (1-3*m_kappa)*dyad(a[0]);
	double E0 = h0.tr() - 1;
	s += h0*(2.*m_k1*E0*exp(m_k2*E0*E0));

	mat3ds h1 = m_kappa*b;
	if (a[1]*a[1] > ONE)
        h1 += (1-3*m_kappa)*dyad(a[1]);
	double E1 = h1.tr() - 1;
	s += h1*(2.*m_k1*E1*exp(m_k2*E1*E1));
    
    return s.dev() / J;
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric tangent
tens4ds FEGasserOgdenHolzapfelUC::DevTangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // determinant of deformation gradient
    double J = pt.m_J;
    
    // Evaluate the distortional deformation gradient
	double Jm13 = pow(J, -1. / 3.);
    mat3d F = pt.m_F*Jm13;
    
    // calculate deviatoric left Cauchy-Green tensor: b = F*Ft
    mat3ds b = pt.LeftCauchyGreen()*(Jm13*Jm13);

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    // Copy the local element basis directions to n
	vec3d n[2];
    n[0].x = Q[0][0]; n[0].y = Q[1][0]; n[0].z = Q[2][0];
    n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1];
    
    // Evaluate the structural direction in the current configuration
    double cg = cos(m_g); double sg = sin(m_g);
	vec3d a[2];
    a[0] = F*(n[0]*cg + n[1]*sg);
    a[1] = F*(n[0]*cg - n[1]*sg);
    
    // Evaluate the ground matrix stress
    mat3ds tau = m_c*b;
    
    // Evaluate the structural tensors in the current configuration
    // and the fiber strains and stress contributions
    mat3ds h0 = m_kappa*b;
    if (a[0]*a[0] > ONE)
        h0 += (1-3*m_kappa)*dyad(a[0]);
	double E0 = h0.tr() - 1;
	double exp0 = exp(m_k2*E0*E0);
	tau += h0*(2.*m_k1*E0*exp0);

	mat3ds h1 = m_kappa*b;
	if (a[1]*a[1] > ONE)
        h1 += (1-3*m_kappa)*dyad(a[1]);
	double E1 = h1.tr() - 1;
	double exp1 = exp(m_k2*E1*E1);
	tau += h1*(2.*m_k1*E1*exp1);

	mat3ds tbar = tau.dev();
    
    // Evaluate the elasticity tensor
    mat3dd I(1);
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    tens4ds C = ((I4 - IxI/3.)*tau.tr()-dyad1s(tbar,I))*(2./3.);
	C += dyad1s(h0.dev())*(4.*m_k1*(1 + 2 * m_k2*E0*E0)*exp0);
	C += dyad1s(h1.dev())*(4.*m_k1*(1 + 2 * m_k2*E1*E1)*exp1);
    
    return C / J;
}

/*
//-----------------------------------------------------------------------------
//! Calculates the deviatoric stress
mat3ds FEGasserOgdenHolzapfelUC::DevStress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    vec3d n[3];			// local element basis directions
    vec3d a[2];			// structural direction in current configuration
    mat3ds h[2];		// structural tensor in current configuration
    double E[2];		// fiber strain
    
    // determinant of deformation gradient
    double J = pt.m_J;
    
    // Evaluate the distortional deformation gradient
    mat3d F = pt.m_F*pow(J,-1./3.);
    
    // calculate deviatoric left Cauchy-Green tensor: b = F*Ft
    mat3ds b = pt.DevLeftCauchyGreen();
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    // Copy the local element basis directions to n
    n[0].x = Q[0][0]; n[0].y = Q[1][0]; n[0].z = Q[2][0];
    n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1];
    
    // Evaluate the structural direction in the current configuration
    double cg = cos(m_g); double sg = sin(m_g);
    a[0] = F*(n[0]*cg + n[1]*sg);
    a[1] = F*(n[0]*cg - n[1]*sg);
    
    // Evaluate the ground matrix stress
    mat3ds s = m_c/J*b;
    
    // Evaluate the structural tensors in the current configuration
    // and the fiber strains and stress contributions
    h[0] = h[1] = m_kappa*b;
    if (a[0]*a[0] > 1) {
        h[0] += (1-3*m_kappa)*dyad(a[0]);
        E[0] = h[0].tr() - 1;
        s += 2./J*m_k1*E[0]*exp(m_k2*E[0]*E[0])*h[0];
    }
    if (a[1]*a[1] > 1) {
        h[1] += (1-3*m_kappa)*dyad(a[1]);
        E[1] = h[1].tr() - 1;
        s += 2./J*m_k1*E[1]*exp(m_k2*E[1]*E[1])*h[1];
    }
    
    return s.dev();
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric tangent
tens4ds FEGasserOgdenHolzapfelUC::DevTangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    vec3d n[3];			// local element basis directions
    vec3d a[2];			// structural direction in current configuration
    mat3ds h[2];		// structural tensor in current configuration
    double E[2];		// fiber strain
    
    // determinant of deformation gradient
    double J = pt.m_J;
    
    // Evaluate the distortional deformation gradient
    mat3d F = pt.m_F*pow(J,-1./3.);
    
    // calculate deviatoric left Cauchy-Green tensor: b = F*Ft
    mat3ds b = pt.DevLeftCauchyGreen();

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    // Copy the local element basis directions to n
    n[0].x = Q[0][0]; n[0].y = Q[1][0]; n[0].z = Q[2][0];
    n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1];
    
    // Evaluate the structural direction in the current configuration
    double cg = cos(m_g); double sg = sin(m_g);
    a[0] = F*(n[0]*cg + n[1]*sg);
    a[1] = F*(n[0]*cg - n[1]*sg);
    
    // Evaluate the ground matrix stress
    mat3ds s = m_c/J*b;
    
    // Evaluate the structural tensors in the current configuration
    // and the fiber strains and stress contributions
    h[0] = h[1] = m_kappa*b;
    if (a[0]*a[0] > 1) {
        h[0] += (1-3*m_kappa)*dyad(a[0]);
        E[0] = h[0].tr() - 1;
        s += 2./J*m_k1*E[0]*exp(m_k2*E[0]*E[0])*h[0];
    }
    if (a[1]*a[1] > 1) {
        h[1] += (1-3*m_kappa)*dyad(a[1]);
        E[1] = h[1].tr() - 1;
        s += 2./J*m_k1*E[1]*exp(m_k2*E[1]*E[1])*h[1];
    }
    
    // Evaluate the elasticity tensor
    mat3dd I(1);
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    tens4ds C = ((I4+IxI/3.)*s.tr()-dyad1s(s,I))*(2./3.);
    if (a[0]*a[0] > 1)
        C += 4./J*m_k1*(1+2*m_k2*E[0]*E[0])*exp(m_k2*E[0]*E[0])*
        (dyad1s(h[0]) - h[0].tr()/3.*(dyad1s(h[0],I) - h[0].tr()/3.*IxI));
    if (a[1]*a[1] > 1)
        C += 4./J*m_k1*(1+2*m_k2*E[1]*E[1])*exp(m_k2*E[1]*E[1])*
        (dyad1s(h[1]) - h[1].tr()/3.*(dyad1s(h[1],I) - h[1].tr()/3.*IxI));
    
    return C;
}
*/

//-----------------------------------------------------------------------------
//! Calculates the deviatoric stress
double FEGasserOgdenHolzapfelUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    vec3d n[3];			// local element basis directions
    vec3d a[2];			// structural direction in current configuration
    mat3ds h[2];		// structural tensor in current configuration
    double E[2];		// fiber strain
    
    // determinant of deformation gradient
    double J = pt.m_J;
    
    // Evaluate the distortional deformation gradient
    mat3d F = pt.m_F*pow(J,-1./3.);
    
    // calculate deviatoric left Cauchy-Green tensor: b = F*Ft
    mat3ds b = pt.DevLeftCauchyGreen();
    double I1 = b.tr();
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    // Copy the local element basis directions to n
    n[0].x = Q[0][0]; n[0].y = Q[1][0]; n[0].z = Q[2][0];
    n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1];
    
    // Evaluate the structural direction in the current configuration
    double cg = cos(m_g); double sg = sin(m_g);
    a[0] = F*(n[0]*cg + n[1]*sg);
    a[1] = F*(n[0]*cg - n[1]*sg);
    
    // Evaluate the ground matrix strain energy density
    double sed = 0.5*m_c*(I1 - 3);
    
    // Evaluate the structural tensors in the current configuration
    // and the fiber strains and stress contributions
    h[0] = h[1] = m_kappa*b;
    if (a[0]*a[0] > 1) {
        h[0] += (1-3*m_kappa)*dyad(a[0]);
        E[0] = h[0].tr() - 1;
        sed += 0.5*m_k1/m_k2*(exp(m_k2*E[0]*E[0])-1);
    }
    if (a[1]*a[1] > 1) {
        h[1] += (1-3*m_kappa)*dyad(a[1]);
        E[1] = h[1].tr() - 1;
        sed += 0.5*m_k1/m_k2*(exp(m_k2*E[1]*E[1])-1);
    }
    
    return sed;
}
