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
#include "FENaturalNeoHookean.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FENaturalNeoHookean, FEElasticMaterial)
    ADD_PARAMETER(m_E, FE_RANGE_GREATER_OR_EQUAL(0.0), "E")->setUnits(UNIT_PRESSURE)->setLongName("Young's modulus");
    ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1, 0.5), "v")->setLongName("Poisson's ratio");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
mat3ds FENaturalNeoHookean::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    double J = pt.m_J;
    
    // get the material parameters
    double E = m_E(mp);
    double v = m_v(mp);
    double k = E/3/(1-2*v);
    double mu = E/2/(1+v);
    
    // evaluate spatial Hencky (logarithmic) strain
    mat3ds h = pt.LeftHencky();
    
    // evaluate amount of dilatation
    double K1 = h.tr();
    
    // evaluate amount of distortion (always positive)
    mat3ds hdev = h.dev();
    double K2 = hdev.norm();
    
    // Identity
    mat3dd I(1);
    
    // Phi
    mat3ds Phi;
    if (K2 != 0) Phi = hdev/K2;
    else Phi.zero();
    
    // calculate stress
    mat3ds s = (I*(k*K1) + Phi*(2*mu*K2))/J;
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FENaturalNeoHookean::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the material parameters
    double E = m_E(mp);
    double v = m_v(mp);
    double k = E/3/(1-2*v);
    double mu = E/2/(1+v);

    // jacobian
    double J = pt.m_J;
    
    // get the left Cauchy-Green tensor
    mat3ds b = pt.LeftCauchyGreen();
    
    // get the eigenvalues and eigenvectors of b
    double lam2[3];    // these are the squares of the eigenvalues of V
    vec3d n[3];
    EigenValues(b, lam2, n, 1e-12);
    
    // get the eigenvalues of V
    double lam[3], lnl[3];
    mat3ds a[3];
    for (int i=0; i<3; ++i) {
        lam[i] = sqrt(lam2[i]);
        lnl[i] = log(lam[i]);
        a[i] = dyad(n[i]);
    }
    
    // evaluate relevant coefficients
    double s[3], ss[3], ddW, ds[3];
    double c1 = 3*k + 4*mu;
    double c2 = 3*k - 2*mu;
    double c3 = 3*k + mu;
    s[0] = (c1*lnl[0] + c2*(lnl[1] + lnl[2]))/(3*J);
    s[1] = (c1*lnl[1] + c2*(lnl[2] + lnl[0]))/(3*J);
    s[2] = (c1*lnl[2] + c2*(lnl[0] + lnl[1]))/(3*J);
    ss[0] = (c1*(1 - 2*lnl[0]) - 2*c2*(lnl[1] + lnl[2]))/(3*J);
    ss[1] = (c1*(1 - 2*lnl[1]) - 2*c2*(lnl[2] + lnl[0]))/(3*J);
    ss[2] = (c1*(1 - 2*lnl[2]) - 2*c2*(lnl[0] + lnl[1]))/(3*J);
    ddW = c2/(3*J);
    if (lam2[1] != lam2[2])
        ds[0] = 2*(lam2[2]*s[1] - lam2[1]*s[2])/(lam2[1] - lam2[2]);
    else
        ds[0] = (6*mu - 4*c3*lnl[1] - 2*c2*lnl[0])/(3*J);
    if (lam2[2] != lam2[0])
        ds[1] = 2*(lam2[0]*s[2] - lam2[2]*s[0])/(lam2[2] - lam2[0]);
    else
        ds[1] = (6*mu - 4*c3*lnl[2] - 2*c2*lnl[1])/(3*J);
    if (lam2[0] != lam2[1])
        ds[2] = 2*(lam2[1]*s[0] - lam2[0]*s[1])/(lam2[0] - lam2[1]);
    else
        ds[2] = (6*mu - 4*c3*lnl[0] - 2*c2*lnl[2])/(3*J);
    
    tens4ds c = dyad1s(a[0])*ss[0] + dyad1s(a[1])*ss[1] + dyad1s(a[2])*ss[2]
    + (dyad1s(a[1], a[2]) + dyad1s(a[2], a[0]) + dyad1s(a[0], a[1]))*ddW
    + dyad4s(a[1], a[2])*ds[0] + dyad4s(a[2], a[0])*ds[1] + dyad4s(a[0], a[1])*ds[2];
    
    return c;
}

//-----------------------------------------------------------------------------
double FENaturalNeoHookean::StrainEnergyDensity(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the material parameters
    double E = m_E(mp);
    double v = m_v(mp);
    double k = E/3/(1-2*v);
    double mu = E/2/(1+v);

    // evaluate spatial Hencky (logarithmic) strain
    mat3ds h = pt.LeftHencky();
    
    // evaluate amount of dilatation
    double K1 = h.tr();
    
    // evaluate amount of distortion (always positive)
    mat3ds hdev = h.dev();
    double K2 = hdev.norm();

    double sed = K1*K1*k/2 + K2*K2*mu;
    
    return sed;
}

//-----------------------------------------------------------------------------
void FENaturalNeoHookean::EigenValues(mat3ds& A, double l[3], vec3d r[3], const double eps)
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

