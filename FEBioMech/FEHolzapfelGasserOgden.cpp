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
#include "FEHolzapfelGasserOgden.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEHolzapfelGasserOgden, FEUncoupledMaterial)
	ADD_PARAMETER(m_c    , FE_RANGE_GREATER_OR_EQUAL(0.0), "c");
	ADD_PARAMETER(m_k1   , FE_RANGE_GREATER_OR_EQUAL(0.0), "k1");
	ADD_PARAMETER(m_k2   , FE_RANGE_GREATER_OR_EQUAL(0.0), "k2");
	ADD_PARAMETER(m_kappa, FE_RANGE_CLOSED(0.0, 1.0/3.0), "kappa");
	ADD_PARAMETER(m_gdeg , "gamma");

    ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEHolzapfelGasserOgden::FEHolzapfelGasserOgden(FEModel* pfem) : FEUncoupledMaterial(pfem) 
{ 
    m_npmodel = 3; 
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric stress
mat3ds FEHolzapfelGasserOgden::DevStress(FEMaterialPoint& mp)
{
    double c = m_c(mp);
    double k1 = m_k1(mp);
    double k2 = m_k2(mp);
    double kappa = m_kappa(mp);
    double g = m_gdeg(mp)*PI/180;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // determinant of deformation gradient
    double J = pt.m_J;

   
    // Evaluate the distortional deformation gradient
	double Jm13 = pow(J, -1. / 3.);
	mat3d F = pt.m_F*Jm13;
    
    // calculate deviatoric left Cauchy-Green tensor: b = F*Ft
    mat3ds b = pt.DevLeftCauchyGreen();
    mat3ds C = pt.DevRightCauchyGreen();
    double I1 = C.tr();
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    // Copy the local element basis directions to n
	vec3d n[2];
    n[0].x = Q[0][0]; n[0].y = Q[1][0]; n[0].z = Q[2][0];
    n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1];
    
    // Evaluate the structural direction in the current configuration
    double cg = cos(g); double sg = sin(g);
	vec3d ar[2],a[2];
    ar[0] = n[0]*cg + n[1]*sg; a[0] = F*ar[0];
    ar[1] = n[0]*cg - n[1]*sg; a[1] = F*ar[1];
    
    // Evaluate the ground matrix stress
    mat3ds s = c*b;
    
    // Evaluate the structural tensors in the current configuration
    // and the fiber strains and stress contributions
    double I40 = ar[0]*(C*ar[0]);
    double E0 = kappa*(I1-3) + (1-3*kappa)*(I40-1);
    if (E0 >= 0) {
        mat3ds h0 = kappa*b + (1-3*kappa)*dyad(a[0]);
        s += h0*(2.*k1*E0*exp(k2*E0*E0));
    }

    double I41 = ar[1]*(C*ar[1]);
    double E1 = kappa*(I1-3) + (1-3*kappa)*(I41-1);
    if (E1 >= 0) {
        mat3ds h1 = kappa*b + (1-3*kappa)*dyad(a[1]);
        s += h1*(2.*k1*E1*exp(k2*E1*E1));
    }

    return s.dev() / J;
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric tangent
tens4ds FEHolzapfelGasserOgden::DevTangent(FEMaterialPoint& mp)
{
    double c = m_c(mp);
    double k1 = m_k1(mp);
    double k2 = m_k2(mp);
    double kappa = m_kappa(mp);
    double g = m_gdeg(mp)*PI/180;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // determinant of deformation gradient
    double J = pt.m_J;
    
    // Evaluate the distortional deformation gradient
	double Jm13 = pow(J, -1. / 3.);
    mat3d F = pt.m_F*Jm13;
    
    // calculate deviatoric left Cauchy-Green tensor: b = F*Ft
    mat3ds b = pt.DevLeftCauchyGreen();
    mat3ds C = pt.DevRightCauchyGreen();
    double I1 = C.tr();

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    // Copy the local element basis directions to n
	vec3d n[2];
    n[0].x = Q[0][0]; n[0].y = Q[1][0]; n[0].z = Q[2][0];
    n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1];
    
    // Evaluate the structural direction in the current configuration
    double cg = cos(g); double sg = sin(g);
    vec3d ar[2],a[2];
    ar[0] = n[0]*cg + n[1]*sg; a[0] = F*ar[0];
    ar[1] = n[0]*cg - n[1]*sg; a[1] = F*ar[1];

    // Evaluate the ground matrix stress
    mat3ds s = c*b;
    
    // Evaluate the structural tensors in the current configuration
    // and the fiber strains and stress contributions
    double I40 = ar[0]*(C*ar[0]);
    double E0 = kappa*(I1-3) + (1-3*kappa)*(I40-1);
    mat3ds h0;
    if (E0 >= 0) {
        h0 = kappa*b + (1-3*kappa)*dyad(a[0]);
        s += h0*(2.*k1*E0*exp(k2*E0*E0));
    }

    double I41 = ar[1]*(C*ar[1]);
    double E1 = kappa*(I1-3) + (1-3*kappa)*(I41-1);
    mat3ds h1;
    if (E1 >= 0) {
        h1 = kappa*b + (1-3*kappa)*dyad(a[1]);
        s += h1*(2.*k1*E1*exp(k2*E1*E1));
    }

	mat3ds sbar = s.dev();
    
    // Evaluate the elasticity tensor
    mat3dd I(1);
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    tens4ds ce = ((I4 - IxI/3.)*s.tr()-dyad1s(sbar,I))*(2./3.);
	if (E0 >= 0) ce += dyad1s(h0.dev())*(4.*k1*(1 + 2 * k2*E0*E0)*exp(k2*E0*E0));
	if (E1 >= 0) ce += dyad1s(h1.dev())*(4.*k1*(1 + 2 * k2*E1*E1)*exp(k2*E1*E1));
    
    return ce / J;
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric stress
double FEHolzapfelGasserOgden::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    double c = m_c(mp);
    double k1 = m_k1(mp);
    double k2 = m_k2(mp);
    double kappa = m_kappa(mp);
    double g = m_gdeg(mp)*PI/180;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    vec3d n[3];			// local element basis directions
    mat3ds h[2];		// structural tensor in current configuration
    double E[2];		// fiber strain
    
    // calculate deviatoric right Cauchy-Green tensor
    mat3ds C = pt.DevRightCauchyGreen();
    double I1 = C.tr();
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    // Copy the local element basis directions to n
    n[0].x = Q[0][0]; n[0].y = Q[1][0]; n[0].z = Q[2][0];
    n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1];
    
    // Evaluate the structural direction in the current configuration
    double cg = cos(g); double sg = sin(g);
    vec3d ar[2];
    ar[0] = n[0]*cg + n[1]*sg;
    ar[1] = n[0]*cg - n[1]*sg;

    // Evaluate the ground matrix strain energy density
    double sed = 0.5*c*(I1 - 3);
    
    // Evaluate the structural tensors in the current configuration
    // and the fiber strains and stress contributions
    double I40 = ar[0]*(C*ar[0]);
    E[0] = kappa*(I1-3) + (1-3*kappa)*(I40-1);
    if (E[0] >= 0)
        sed += 0.5*k1/k2*(exp(k2*E[0]*E[0])-1);
    double I41 = ar[1]*(C*ar[1]);
    E[1] = kappa*(I1-3) + (1-3*kappa)*(I41-1);
    if (E[1] >= 0)
        sed += 0.5*k1/k2*(exp(k2*E[1]*E[1])-1);
    
    return sed;
}
