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
#include "FECarreauYasudaViscousSolid.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FECarreauYasudaViscousSolid, FEElasticMaterial)
    ADD_PARAMETER(m_mu0, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu0");
    ADD_PARAMETER(m_mui, FE_RANGE_GREATER_OR_EQUAL(0.0), "mui");
    ADD_PARAMETER(m_lam, FE_RANGE_GREATER_OR_EQUAL(0.0), "lambda");
    ADD_PARAMETER(m_n  , FE_RANGE_GREATER_OR_EQUAL(0.0), "n");
    ADD_PARAMETER(m_a  , FE_RANGE_GREATER_OR_EQUAL(2.0), "a");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
mat3ds FECarreauYasudaViscousSolid::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    mat3ds D = pt.RateOfDeformation();
//    double d[3];
//    D.eigen(d);
    
    double mu0 = m_mu0(mp);
    double mui = m_mui(mp);
    double lam = m_lam(mp);
    double n = m_n(mp);
    double a = m_a(mp);
    double gdot = sqrt(2*(D.sqr()).tr());
//    double gdot = max(fabs(d[2]-d[1]),max(fabs(d[0]-d[2]),fabs(d[1]-d[0])));
    double lamga = pow(lam*gdot,a);
    double mu = mui + (mu0 - mui)*pow(1+lamga, (n-1)/a);
    mat3ds s = D*(2*mu);

    return s;
}

//-----------------------------------------------------------------------------
tens4ds FECarreauYasudaViscousSolid::Tangent(FEMaterialPoint& mp)
{
    const FETimeInfo& tp = GetTimeInfo();
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    tens4ds Cv;

    if (tp.timeIncrement > 0) {
        mat3ds D = pt.RateOfDeformation();

        double mu0 = m_mu0(mp);
        double mui = m_mui(mp);
        double lam = m_lam(mp);
        double n = m_n(mp);
        double a = m_a(mp);
        double gdot = sqrt(2*(D.sqr()).tr());
        double lamga = pow(lam*gdot,a);
        double mu = mui + (mu0 - mui)*pow(1+lamga, (n-1)/a);
        double dmu = 2*(mu0 - mui)*(n-1)*pow(lam,a)*pow(gdot,a-2)*pow(1+lamga, (n-a-1)/a);

        mat3dd I(1);
        double tmp = tp.alphaf*tp.gamma/(tp.beta*tp.timeIncrement);
        Cv = (dyad1s(D)*(2*dmu) + dyad4s(I)*(2*mu))*tmp;
    }
    else Cv.zero();

    return Cv;
}

//-----------------------------------------------------------------------------
double FECarreauYasudaViscousSolid::StrainEnergyDensity(FEMaterialPoint& mp)
{
    return 0;
}

