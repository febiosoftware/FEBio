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
#include "FEDamageCDF.h"
#include "FEDamageCriterion.h"
#include "FEDamageMaterialPoint.h"
#include <FECore/log.h>
#include <FECore/gamma.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDamageCDF, FEMaterialProperty)
	ADD_PARAMETER(m_Dmax , FE_RANGE_CLOSED(0.0, 1.0), "Dmax");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
double FEDamageCDF::Damage(FEMaterialPoint& mp) {
    
    // get the damage material point data
    FEDamageMaterialPoint& dp = *mp.ExtractData<FEDamageMaterialPoint>();
    
    // get the damage criterion (assuming dp.m_Etrial is up-to-date)
    double Es = max(dp.m_Etrial, dp.m_Emax);

    dp.m_D = cdf(mp,Es)*m_Dmax;
    
    return dp.m_D;
}


///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDamageCDFSimo, FEDamageCDF)
    ADD_PARAMETER(m_alpha, FE_RANGE_GREATER(0.0), "a");
    ADD_PARAMETER(m_beta , FE_RANGE_CLOSED(0.0, 1.0), "b");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFSimo::FEDamageCDFSimo(FEModel* pfem) : FEDamageCDF(pfem)
{
}

//-----------------------------------------------------------------------------
// Simo damage cumulative distribution function
// Simo, CMAME 60 (1987), 153-173
// Simo damage cumulative distribution function
double FEDamageCDFSimo::cdf(FEMaterialPoint& mp, const double X)
{
    double alpha = m_alpha(mp);
    double beta = m_beta(mp);
    if (alpha == 0) {
        return 0;
    }
    
    double cdf = 0;
    
    // this CDF only admits positive values
    if (X >= 0) {
        if (X > 1e-12) cdf = 1 - beta - (1.0 - beta)*(1.0 - exp(-X/alpha))*alpha/X;
        else cdf = 0.5*(1.0 - beta)/alpha*X;
    }
    
    return cdf;
}

// Simo damage probability density function
double FEDamageCDFSimo::pdf(FEMaterialPoint& mp, const double X)
{
    double alpha = m_alpha(mp);
    double beta = m_beta(mp);
    if (alpha == 0) {
        return 0;
    }
    
    double pdf = 0;
    
    // this CDF only admits positive values
    if (X >= 0) {
        if (X > 1e-12) pdf = (1.0 - beta)*(alpha - (alpha + X)*exp(-X/alpha))/(X*X);
        else pdf = (1.0 - beta)/alpha*(0.5 - X/3/alpha);
    }
    
    return pdf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDamageCDFLogNormal, FEDamageCDF)
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER(0.0), "mu");
    ADD_PARAMETER(m_sigma, FE_RANGE_GREATER(0.0), "sigma");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFLogNormal::FEDamageCDFLogNormal(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_mu = 1;
    m_sigma = 1;
    m_Dmax = 1;
}

//-----------------------------------------------------------------------------
// Lognormal damage cumulative distribution function
double FEDamageCDFLogNormal::cdf(FEMaterialPoint& mp, const double X)
{
    double cdf = 0;
    double mu = m_mu(mp);
    double sigma = m_sigma(mp);
    // this CDF only admits positive values
    if (X >= 0)
        cdf = 0.5*erfc(-log(X/mu)/sigma/sqrt(2.));
    
    return cdf;
}

// Lognormal damage probability density function
double FEDamageCDFLogNormal::pdf(FEMaterialPoint& mp, const double X)
{
    double pdf = 0;
    double mu = m_mu(mp);
    double sigma = m_sigma(mp);

    // this CDF only admits positive values
    if (X > 1e-12) pdf = exp(-pow(log(X/mu)/sigma,2)/2)/(sqrt(2*M_PI)*X*sigma);
    
    return pdf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDamageCDFWeibull, FEDamageCDF)
    ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
    ADD_PARAMETER(m_ploc, "ploc");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFWeibull::FEDamageCDFWeibull(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_alpha = 1;
    m_mu = 1;
    m_ploc = 0;
}

//-----------------------------------------------------------------------------
// Weibull damage cumulative distribution function
double FEDamageCDFWeibull::cdf(FEMaterialPoint& mp, const double X)
{
    double cdf = 0;
    double alpha = m_alpha(mp);
    double mu = m_mu(mp);
    double ploc = m_ploc(mp);
    
    // this CDF only admits positive values
    if (X > ploc)
        cdf = 1 - exp(-pow((X-ploc)/mu,alpha));
    
    return cdf;
}

// Weibull damage probability density function
double FEDamageCDFWeibull::pdf(FEMaterialPoint& mp, const double X)
{
    double pdf = 0;
    double alpha = m_alpha(mp);
    double mu = m_mu(mp);
    double ploc = m_ploc(mp);

    // this CDF only admits positive values
    if ((alpha > 1) && (X > ploc))
        pdf = exp(-pow((X-ploc)/mu,alpha))*alpha*pow(X-ploc, alpha-1)/pow(mu, alpha);
    else if ((alpha == 1) && (X > ploc))
        pdf = exp(-(X-ploc)/mu)/mu;
    
    return pdf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDamageCDFStep, FEDamageCDF)
    ADD_PARAMETER(m_mu  , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFStep::FEDamageCDFStep(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_mu = 1;
}

//-----------------------------------------------------------------------------
// Step cumulative distribution function (sudden fracture)
// Step damage cumulative distribution function
double FEDamageCDFStep::cdf(FEMaterialPoint& mp, const double X)
{
    double cdf = 0;
    double mu = m_mu(mp);
    
    // this CDF only admits positive values
    if (X > mu)
        cdf = 1.0;
    
    return cdf;
}

// Step damage probability density function
double FEDamageCDFStep::pdf(FEMaterialPoint& mp, const double X)
{
    double pdf = 0;
    double mu = m_mu(mp);

    // this PDF only admits positive values
    if (X == mu) pdf = 1.0;
    
    return pdf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDamageCDFPQP, FEDamageCDF)
    ADD_PARAMETER(m_mumin, FE_RANGE_GREATER_OR_EQUAL(0.0), "mumin");
    ADD_PARAMETER(m_mumax, "mumax");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFPQP::FEDamageCDFPQP(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_mumin = 0;
    m_mumax = 1;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEDamageCDFPQP::Validate()
{
	return FEDamageCDF::Validate();
}

//-----------------------------------------------------------------------------
// Piecewise S-shaped quintic polynomial damage cumulative distribution function
double FEDamageCDFPQP::cdf(FEMaterialPoint& mp, const double X)
{
    double cdf = 0;
    double mumin = m_mumin(mp);
    double mumax = m_mumax(mp);

    if (X <= mumin) cdf = 0;
    else if (X >= mumax) cdf = 1;
    else
    {
        double x = (X - mumin)/(mumax - mumin);
        cdf = pow(x,3)*(10 - 15*x + 6*x*x);
    }

    return cdf;
}

// Piecewise S-shaped quintic polynomial damage probability density function
double FEDamageCDFPQP::pdf(FEMaterialPoint& mp, const double X)
{
    double pdf = 0;
    double mumin = m_mumin(mp);
    double mumax = m_mumax(mp);

    if (X <= mumin) pdf = 0;
    else if (X >= mumax) pdf = 0;
    else
    {
        double x = (X - mumin)/(mumax - mumin);
        pdf = pow(x*(x-1),2)*30;
    }

    return pdf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDamageCDFGamma, FEDamageCDF)
    ADD_PARAMETER(m_alpha, FE_RANGE_GREATER(0)           , "alpha");
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFGamma::FEDamageCDFGamma(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_alpha = 2;
    m_mu = 4;
}

//-----------------------------------------------------------------------------
// Gamma damage cumulative distribution function
double FEDamageCDFGamma::cdf(FEMaterialPoint& mp, const double X)
{
    double cdf = 0;
    double alpha = m_alpha(mp);
    double mu = m_mu(mp);
    
    // this CDF only admits positive values
    if (X > 0)
        cdf = gamma_inc_P(alpha,X/mu);
    
    return cdf;
}

// Gamma damage probability density function
double FEDamageCDFGamma::pdf(FEMaterialPoint& mp, const double X)
{
    double pdf = 0;
    double alpha = m_alpha(mp);
    double mu = m_mu(mp);

    // this CDF only admits positive values
    if (X > 0)
        pdf = pow(X/mu, alpha-1)*exp(-X/mu)/mu*gammainv(alpha);
    
    return pdf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDamageCDFUser, FEDamageCDF)
    ADD_PROPERTY(m_cdf, "cdf");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFUser::FEDamageCDFUser(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_cdf = nullptr;
}

//-----------------------------------------------------------------------------
// User-defined loadcurve for damage cumulative distribution function
double FEDamageCDFUser::cdf(FEMaterialPoint& mp, const double X)
{
    return m_cdf->value(X);
}

// Derivative of user-defined loadcurve damage probability density function
double FEDamageCDFUser::pdf(FEMaterialPoint& mp, const double X)
{
    return m_cdf->derive(X);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDamageCDFPower, FEDamageCDF)
    ADD_PARAMETER(m_alpha, "alpha" )->setLongName("power exponent");
    ADD_PARAMETER(m_mu0  , "mu0"   )->setLongName("constant mu0");
    ADD_PARAMETER(m_mu1  , "mu1"   )->setLongName("power coefficient");
    ADD_PARAMETER(m_s    , FE_RANGE_GREATER(0.0)         , "scale" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFPower::FEDamageCDFPower(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_alpha = 2;
    m_mu0 = 1;
    m_mu1 = 0;
    m_s = 1;
}

//-----------------------------------------------------------------------------
// Power cumulative distribution function
double FEDamageCDFPower::cdf(FEMaterialPoint& mp, const double X)
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

// Power probability density function
double FEDamageCDFPower::pdf(FEMaterialPoint& mp, const double X)
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
BEGIN_FECORE_CLASS(FEDamageCDFExp, FEDamageCDF)
    ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(1.0), "alpha" );
    ADD_PARAMETER(m_mu0  , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu0"   );
    ADD_PARAMETER(m_mu1  , "mu1"   );
    ADD_PARAMETER(m_s    , FE_RANGE_GREATER(0.0)         , "scale" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFExp::FEDamageCDFExp(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_alpha = 2;
    m_mu0 = 1;
    m_mu1 = 0;
    m_s = 1;
}

//-----------------------------------------------------------------------------
// Exponential cumulative distribution function
double FEDamageCDFExp::cdf(FEMaterialPoint& mp, const double X)
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
double FEDamageCDFExp::pdf(FEMaterialPoint& mp, const double X)
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
BEGIN_FECORE_CLASS(FEDamageCDFPoly2, FEDamageCDF)
    ADD_PARAMETER(m_mu0  , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu0"   );
    ADD_PARAMETER(m_mu1  , "mu1"   );
    ADD_PARAMETER(m_mu2  , "mu2"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFPoly2::FEDamageCDFPoly2(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_mu0 = 1;
    m_mu1 = 0;
    m_mu2 = 0;
}

//-----------------------------------------------------------------------------
// Poly2 cumulative distribution function
double FEDamageCDFPoly2::cdf(FEMaterialPoint& mp, const double X)
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
double FEDamageCDFPoly2::pdf(FEMaterialPoint& mp, const double X)
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

