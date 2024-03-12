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
    ADD_PROPERTY(m_brf, "function");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRecruitmentUser::FEBondRecruitmentUser(FEModel* pfem) : FEBondRecruitment(pfem)
{
    m_brf = nullptr;
}

//-----------------------------------------------------------------------------
// User-defined loadcurve for bond recruitment function
double FEBondRecruitmentUser::brf(FEMaterialPoint& mp, const double X)
{
    return m_brf->value(X);
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
double FEBondRecruitmentPower::brf(FEMaterialPoint& mp, const double X)
{
    double brf = 0;
    double alpha = m_alpha(mp);
    double mu0 = m_mu0(mp);
    double mu1 = m_mu1(mp);
    double s = m_s(mp);

    // this brf only admits positive values
    if (X > 0)
        brf = mu0 + mu1*pow(X/s,alpha);
    
    return brf;
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
double FEBondRecruitmentExp::brf(FEMaterialPoint& mp, const double X)
{
    double brf = 0;
    double alpha = m_alpha(mp);
    double mu0 = m_mu0(mp);
    double mu1 = m_mu1(mp);
    double s = m_s(mp);

    // this brf only admits positive values
    if (X > 0)
        brf = mu0*exp(mu1*pow(X/s,alpha));
    
    return brf;
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
double FEBondRecruitmentPoly::brf(FEMaterialPoint& mp, const double X)
{
    double brf = 0;
    double mu0 = m_mu0(mp);
    double mu1 = m_mu1(mp);
    double mu2 = m_mu2(mp);

    // this brf only admits positive arguments
    if (X > 0)
        brf = mu0 + mu1*X + mu2*X*X;
    
    return brf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEBondRecruitmentLogNormal, FEBondRecruitment)
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER(0.0), "mu");
    ADD_PARAMETER(m_sigma, FE_RANGE_GREATER(0.0), "sigma");
    ADD_PARAMETER(m_max  , FE_RANGE_GREATER_OR_EQUAL(0.0), "max")->setLongName("max recruitment increase");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRecruitmentLogNormal::FEBondRecruitmentLogNormal(FEModel* pfem) : FEBondRecruitment(pfem)
{
    m_mu = 1;
    m_sigma = 1;
    m_max = 0;
}

//-----------------------------------------------------------------------------
// Lognormal damage cumulative distribution function
double FEBondRecruitmentLogNormal::brf(FEMaterialPoint& mp, const double X)
{
    double brf = 0;
    double mu = m_mu(mp);
    double sigma = m_sigma(mp);
    double maxr = m_max(mp);
    // this brf only admits positive values
    if (X >= 0)
        brf = 1.0 + maxr*(0.5*erfc(-log(X/mu)/sigma/sqrt(2.)));
    
    return brf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEBondRecruitmentWeibull, FEBondRecruitment)
    ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
    ADD_PARAMETER(m_ploc, "ploc");
    ADD_PARAMETER(m_max  , FE_RANGE_GREATER_OR_EQUAL(0.0), "max")->setLongName("max recruitment increase");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRecruitmentWeibull::FEBondRecruitmentWeibull(FEModel* pfem) : FEBondRecruitment(pfem)
{
    m_alpha = 1;
    m_mu = 1;
    m_ploc = 0;
    m_max = 0;
}

//-----------------------------------------------------------------------------
// Weibull damage cumulative distribution function
double FEBondRecruitmentWeibull::brf(FEMaterialPoint& mp, const double X)
{
    double brf = 0;
    double alpha = m_alpha(mp);
    double mu = m_mu(mp);
    double ploc = m_ploc(mp);
    double maxr = m_max(mp);

    // this brf only admits positive values
    if (X > ploc)
        brf = 1 + maxr*(1 - exp(-pow((X-ploc)/mu,alpha)));
    
    return brf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEBondRecruitmentPQP, FEBondRecruitment)
    ADD_PARAMETER(m_mumin, FE_RANGE_GREATER_OR_EQUAL(0.0), "mumin");
    ADD_PARAMETER(m_mumax, "mumax");
    ADD_PARAMETER(m_max  , FE_RANGE_GREATER_OR_EQUAL(0.0), "max")->setLongName("max recruitment increase");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRecruitmentPQP::FEBondRecruitmentPQP(FEModel* pfem) : FEBondRecruitment(pfem)
{
    m_mumin = 0;
    m_mumax = 1;
    m_max = 0;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEBondRecruitmentPQP::Validate()
{
    return FEBondRecruitment::Validate();
}

//-----------------------------------------------------------------------------
// Piecewise S-shaped quintic polynomial damage bond recruitment function
double FEBondRecruitmentPQP::brf(FEMaterialPoint& mp, const double X)
{
    double brf = 0;
    double mumin = m_mumin(mp);
    double mumax = m_mumax(mp);
    double maxr = m_max(mp);

    if (X <= mumin) brf = 1;
    else if (X >= mumax) brf = 1+maxr;
    else
    {
        double x = (X - mumin)/(mumax - mumin);
        brf = 1 + maxr*(pow(x,3)*(10 - 15*x + 6*x*x));
    }
    
    return brf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEBondRecruitmentGamma, FEBondRecruitment)
    ADD_PARAMETER(m_alpha, FE_RANGE_GREATER(0), "alpha");
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"   );
    ADD_PARAMETER(m_max  , FE_RANGE_GREATER_OR_EQUAL(0.0), "max")->setLongName("max recruitment increase");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRecruitmentGamma::FEBondRecruitmentGamma(FEModel* pfem) : FEBondRecruitment(pfem)
{
    m_alpha = 2;
    m_mu = 4;
    m_max = 0;
}

//-----------------------------------------------------------------------------
// Gamma damage cumulative distribution function
double FEBondRecruitmentGamma::brf(FEMaterialPoint& mp, const double X)
{
    double brf = 1;
    double alpha = m_alpha(mp);
    double mu = m_mu(mp);
    double maxr = m_max(mp);

    // this CDF only admits positive values
    double scale = gamma_inc_Q(alpha,0);
    if (X > 0)
        brf = 1 + maxr*gamma_inc_P(alpha,X/mu)/scale;
    
    return brf;
}
