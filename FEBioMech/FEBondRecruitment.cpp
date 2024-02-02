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
#include "FEBondRecruitment.h"
#include <FECore/log.h>
#include <FECore/gamma.h>
#include <math.h>

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEBondRecruitmentUser, FEBondRecruitment)
    ADD_PROPERTY(m_cdf, "cdf");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRecruitmentUser::FEBondRecruitmentUser(FEModel* pfem) : FEBondRecruitment(pfem)
{
    m_cdf = nullptr;
}

//-----------------------------------------------------------------------------
// User-defined loadcurve for bond recruitment function
double FEBondRecruitmentUser::cdf(FEMaterialPoint& mp, const double X)
{
    return m_cdf->value(X);
}

// Derivative of user-defined loadcurve damage probability density function
double FEBondRecruitmentUser::pdf(FEMaterialPoint& mp, const double X)
{
    return m_cdf->derive(X);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEBondRecruitmentPower, FEBondRecruitment)
    ADD_PARAMETER(m_alpha, "alpha" )->setLongName("power exponent");
    ADD_PARAMETER(m_mu0  , "mu0"   )->setLongName("constant mu0");
    ADD_PARAMETER(m_mu1  , "mu1"   )->setLongName("power coefficient");
    ADD_PARAMETER(m_s    , FE_RANGE_GREATER(0.0)         , "scale" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRecruitmentPower::FEBondRecruitmentPower(FEModel* pfem) : FEBondRecruitment(pfem)
{
    m_alpha = 2;
    m_mu0 = 1;
    m_mu1 = 0;
    m_s = 1;
}

//-----------------------------------------------------------------------------
// Power bond recruitment function
double FEBondRecruitmentPower::cdf(FEMaterialPoint& mp, const double X)
{
    double cdf = 0;
    double alpha = m_alpha(mp);
    double mu0 = m_mu0(mp);
    double mu1 = m_mu1(mp);
    double s = m_s(mp);

    // this CDF only admits positive values
    if (X > 0)
        cdf = mu0 + mu1*pow(X/s,alpha);
    
    return cdf;
}

// Power bond recruitment function
double FEBondRecruitmentPower::pdf(FEMaterialPoint& mp, const double X)
{
    double pdf = 0;
    double alpha = m_alpha(mp);
    double mu0 = m_mu0(mp);
    double mu1 = m_mu1(mp);
    double s = m_s(mp);

    // this CDF only admits positive values
    if (X > 0)
        pdf = mu1*alpha/s*pow(X/s,alpha-1);
    
    return pdf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEBondRecruitmentExp, FEBondRecruitment)
    ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha" );
    ADD_PARAMETER(m_mu0  , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu0"   );
    ADD_PARAMETER(m_mu1  , "mu1"   );
    ADD_PARAMETER(m_s    , FE_RANGE_GREATER(0.0)         , "scale" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRecruitmentExp::FEBondRecruitmentExp(FEModel* pfem) : FEBondRecruitment(pfem)
{
    m_alpha = 1;
    m_mu0 = 1;
    m_mu1 = 0;
    m_s = 1;
}

//-----------------------------------------------------------------------------
// Exponential bond recruitment function
double FEBondRecruitmentExp::cdf(FEMaterialPoint& mp, const double X)
{
    double cdf = 0;
    double alpha = m_alpha(mp);
    double mu0 = m_mu0(mp);
    double mu1 = m_mu1(mp);
    double s = m_s(mp);

    // this CDF only admits positive values
    if (X > 0)
        cdf = mu0*exp(mu1*pow(X/s,alpha));
    
    return cdf;
}

// Exponential probability density function
double FEBondRecruitmentExp::pdf(FEMaterialPoint& mp, const double X)
{
    double pdf = 0;
    double alpha = m_alpha(mp);
    double mu0 = m_mu0(mp);
    double mu1 = m_mu1(mp);
    double s = m_s(mp);

    // this CDF only admits positive values
    if (X > 0)
        pdf = mu0*mu1*alpha/s*pow(X/s,alpha-1)*exp(mu1*pow(X/s,alpha));
    
    return pdf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEBondRecruitmentPoly, FEBondRecruitment)
    ADD_PARAMETER(m_mu0  , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu0"   );
    ADD_PARAMETER(m_mu1  , "mu1"   );
    ADD_PARAMETER(m_mu2  , "mu2"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRecruitmentPoly::FEBondRecruitmentPoly(FEModel* pfem) : FEBondRecruitment(pfem)
{
    m_mu0 = 1;
    m_mu1 = 0;
    m_mu2 = 0;
}

//-----------------------------------------------------------------------------
// Poly bond recruitment function
double FEBondRecruitmentPoly::cdf(FEMaterialPoint& mp, const double X)
{
    double cdf = 0;
    double mu0 = m_mu0(mp);
    double mu1 = m_mu1(mp);
    double mu2 = m_mu2(mp);

    // this CDF only admits positive arguments
    if (X > 0)
        cdf = mu0 + mu1*X + mu2*X*X;
    
    return cdf;
}

// Poly2 probability density function
double FEBondRecruitmentPoly::pdf(FEMaterialPoint& mp, const double X)
{
    double pdf = 0;
    double mu0 = m_mu0(mp);
    double mu1 = m_mu1(mp);
    double mu2 = m_mu2(mp);

    // this PDF only admits positive arguments
    if (X > 0)
        pdf = mu1 + 2*mu2*X;
    
    return pdf;
}

