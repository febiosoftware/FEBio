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
#include "FEOgdenMaterial.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEOgdenMaterial, FEUncoupledMaterial);
	ADD_PARAMETER(m_c[0], "c1")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c[1], "c2")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c[2], "c3")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c[3], "c4")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c[4], "c5")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c[5], "c6")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_m[0], FE_RANGE_NOT_EQUAL(0.0), "m1");
	ADD_PARAMETER(m_m[1], FE_RANGE_NOT_EQUAL(0.0), "m2");
	ADD_PARAMETER(m_m[2], FE_RANGE_NOT_EQUAL(0.0), "m3");
	ADD_PARAMETER(m_m[3], FE_RANGE_NOT_EQUAL(0.0), "m4");
	ADD_PARAMETER(m_m[4], FE_RANGE_NOT_EQUAL(0.0), "m5");
	ADD_PARAMETER(m_m[5], FE_RANGE_NOT_EQUAL(0.0), "m6");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEOgdenMaterial::FEOgdenMaterial(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	for (int i=0; i<MAX_TERMS; ++i)
	{
		m_c[i] = 0;
		m_m[i] = 1;
	}

	m_eps = 1e-12;
}

//-----------------------------------------------------------------------------
void FEOgdenMaterial::EigenValues(mat3ds& A, double l[3], vec3d r[3], const double eps)
{
    A.eigen(l, r);
    
    // correct for numerical inaccuracy
    double d01 = fabs(l[0] - l[1]);
    double d12 = fabs(l[1] - l[2]);
    double d02 = fabs(l[0] - l[2]);
    
    if (d01 < eps) l[1] = l[0]; //= 0.5*(l[0]+l[1]);
    if (d02 < eps) l[2] = l[0]; //= 0.5*(l[0]+l[2]);
    if (d12 < eps) l[2] = l[1]; //= 0.5*(l[1]+l[2]);
    
}

//-----------------------------------------------------------------------------
//! Calculates the Cauchy stress
mat3ds FEOgdenMaterial::DevStress(FEMaterialPoint &mp)
{
    // extract elastic material data
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // jacobian
    double J = pt.m_J;
    
    // get the left Cauchy-Green tensor
    mat3ds b = pt.DevLeftCauchyGreen();
    
    // get the eigenvalues and eigenvectors of b
    double lam2[3];	// these are the squares of the eigenvalues of V
    vec3d ev[3];
    EigenValues(b, lam2, ev, m_eps);
    
    // get the eigenvalues of V
    double lam[3];
    lam[0] = sqrt(lam2[0]);
    lam[1] = sqrt(lam2[1]);
    lam[2] = sqrt(lam2[2]);

    double c[MAX_TERMS] = { 0 };
    for (int i = 0; i < MAX_TERMS; ++i) c[i] = m_c[i](mp);
    
    // stress
    mat3ds s;
    s.zero();
    double T;
    for (int i=0; i<3; ++i) {
        T = 0;
        for (int j=0; j<MAX_TERMS; ++j)
            T += c[j]/m_m[j]*(pow(lam[i], m_m[j]) - 1)/J;
        s += dyad(ev[i])*T;
    }
    
    return s.dev();
}

//-----------------------------------------------------------------------------
//! Calculates the spatial tangent
tens4ds FEOgdenMaterial::DevTangent(FEMaterialPoint& mp)
{
    int i,j,k;
    
    // extract elastic material data
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // jacobian
    double J = pt.m_J;
    
    // get the left Cauchy-Green tensor
    mat3ds b = pt.DevLeftCauchyGreen();
    
    // get the eigenvalues and eigenvectors of b
    double lam2[3];	// these are the squares of the eigenvalues of V
    vec3d ev[3];
    EigenValues(b, lam2, ev, m_eps);
    
    // get the eigenvalues of V
    double lam[3];
    mat3ds N[3];
    for (i=0; i<3; ++i) {
        lam[i] = sqrt(lam2[i]);
        N[i] = dyad(ev[i]);
    }
    
    // calculate the powers of eigenvalues
    double lamp[3][MAX_TERMS];
    for (j=0; j<MAX_TERMS; ++j)
    {
        lamp[0][j] = pow(lam[0], m_m[j]);
        lamp[1][j] = pow(lam[1], m_m[j]);
        lamp[2][j] = pow(lam[2], m_m[j]);
    }

    // evaluate coefficients at material points
    double ci[MAX_TERMS] = { 0 };
    for (int i = 0; i < MAX_TERMS; ++i) ci[i] = m_c[i](mp);
    
    // principal stresses
    mat3ds s;
    s.zero();
    double T[3];
    for (i=0; i<3; ++i) {
        T[i] = 0;
        for (j=0; j<MAX_TERMS; ++j)
            T[i] += ci[j]/m_m[j]*(lamp[i][j] - 1)/J;
        s += N[i]*T[i];
    }
    
    // coefficients appearing in elasticity tensor
    double D[3][3],E[3][3];
    for (j=0; j<3; ++j) {
        D[j][j] = 0;
        for (k=0; k<MAX_TERMS; ++k)
            D[j][j] += ci[k]/m_m[k]*((m_m[k]-2)*lamp[j][k]+2)/J;
        for (i=j+1; i<3; ++i) {
            D[i][j] = 0;
            if (lam2[j] != lam2[i])
                E[i][j] = 2*(lam2[j]*T[i] - lam2[i]*T[j])/(lam2[i]-lam2[j]);
            else {
                E[i][j] = 0;
                for (k=0; k<MAX_TERMS; ++k)
                    E[i][j] += ci[k]/m_m[k]*((m_m[k]-2)*lamp[j][k]+2)/J;
            }
        }
    }
    
    // spatial elasticity tensor
    tens4ds c(0.0);
    mat3dd I(1.0);
    for (j=0; j<3; ++j) {
        c += dyad1s(N[j])*D[j][j];
        for (i=j+1; i<3; ++i) {
            c += dyad1s(N[i],N[j])*D[i][j];
            c += dyad4s(N[i],N[j])*E[i][j];
        }
    }
    
    // This is the final value of the elasticity tensor
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    
    c += - 1./3.*(ddots(c,IxI) - IxI*(c.tr()/3.))
    + 2./3.*((I4-IxI/3.)*s.tr()-dyad1s(s.dev(),I));
    
    return c;
}

//-----------------------------------------------------------------------------
double FEOgdenMaterial::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    // extract elastic material data
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the left Cauchy-Green tensor
    mat3ds b = pt.DevLeftCauchyGreen();
    
    // get the eigenvalues and eigenvectors of b
    double lam2[3];	// these are the squares of the eigenvalues of V
    vec3d ev[3];
    EigenValues(b, lam2, ev, m_eps);
    
    // get the eigenvalues of V
    double lam[3];
    lam[0] = sqrt(lam2[0]);
    lam[1] = sqrt(lam2[1]);
    lam[2] = sqrt(lam2[2]);
    
    // evaluate coefficients at material points
    double ci[MAX_TERMS] = { 0 };
    for (int i = 0; i < MAX_TERMS; ++i) ci[i] = m_c[i](mp);

    // strain energy density
    double sed = 0;
    for (int j=0; j<MAX_TERMS; ++j)
        sed += ci[j]/(m_m[j]*m_m[j])*
        (pow(lam[0], m_m[j]) + pow(lam[1], m_m[j]) + pow(lam[2], m_m[j]) - 3);
    
    return sed;
}
