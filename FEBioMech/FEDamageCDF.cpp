/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#define _USE_MATH_DEFINES
#include <math.h>
#ifdef HAVE_GSL
#include "gsl/gsl_sf_gamma.h"
#endif

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDamageCDF, FEMaterial)
	ADD_PARAMETER(m_Dmax , FE_RANGE_CLOSED(0.0, 1.0), "Dmax");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
double FEDamageCDF::Damage(FEMaterialPoint& mp) {
    
    // get the damage material point data
    FEDamageMaterialPoint& dp = *mp.ExtractData<FEDamageMaterialPoint>();
    
    // get the damage criterion (assuming dp.m_Etrial is up-to-date)
    double Es = max(dp.m_Etrial, dp.m_Emax);

    dp.m_D = cdf(Es)*m_Dmax;
    
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
double FEDamageCDFSimo::cdf(const double X)
{
    if (m_alpha == 0) {
        return 0;
    }
    
    double cdf = 0;
    
    // this CDF only admits positive values
    if (X >= 0) {
        if (X > 1e-12) cdf = 1 - m_beta - (1.0 - m_beta)*(1.0 - exp(-X/m_alpha))*m_alpha/X;
        else cdf = 0.5*(1.0 - m_beta)/m_alpha*X;
    }
    
    return cdf;
}

// Simo damage probability density function
double FEDamageCDFSimo::pdf(const double X)
{
    if (m_alpha == 0) {
        return 0;
    }
    
    double pdf = 0;
    
    // this CDF only admits positive values
    if (X >= 0) {
        if (X > 1e-12) pdf = (1.0 - m_beta)*(m_alpha - (m_alpha + X)*exp(-X/m_alpha))/(X*X);
        else pdf = (1.0 - m_beta)/m_alpha*(0.5 - X/3/m_alpha);
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
double FEDamageCDFLogNormal::cdf(const double X)
{
    double cdf = 0;
    
    // this CDF only admits positive values
    if (X >= 0)
        cdf = 0.5*erfc(-log(X/m_mu)/m_sigma/sqrt(2.));
    
    return cdf;
}

// Lognormal damage probability density function
double FEDamageCDFLogNormal::pdf(const double X)
{
    double pdf = 0;
    
    // this CDF only admits positive values
    if (X > 1e-12) pdf = exp(-pow(log(X/m_mu)/m_sigma,2)/2)/(sqrt(2*M_PI)*X*m_sigma);
    
    return pdf;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDamageCDFWeibull, FEDamageCDF)
    ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFWeibull::FEDamageCDFWeibull(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_alpha = m_mu;
}

//-----------------------------------------------------------------------------
// Weibull damage cumulative distribution function
double FEDamageCDFWeibull::cdf(const double X)
{
    double cdf = 0;
    
    // this CDF only admits positive values
    if (X > 0)
        cdf = 1 - exp(-pow(X/m_mu,m_alpha));
    
    return cdf;
}

// Weibull damage probability density function
double FEDamageCDFWeibull::pdf(const double X)
{
    double pdf = 0;
    
    // this CDF only admits positive values
    if ((m_alpha > 1) && (X > 0))
        pdf = exp(-pow(X/m_mu,m_alpha))*m_alpha*pow(X, m_alpha-1)/pow(m_mu, m_alpha);
    else if (m_alpha == 1)
        pdf = exp(-X/m_mu)/m_mu;
    
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
double FEDamageCDFStep::cdf(const double X)
{
    double cdf = 0;
    
    // this CDF only admits positive values
    if (X > m_mu)
        cdf = 1.0;
    
    return cdf;
}

// Step damage probability density function
double FEDamageCDFStep::pdf(const double X)
{
    double pdf = 0;
    
    // this PDF only admits positive values
    if (X == m_mu) pdf = 1.0;
    
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
	if (m_mumax <= m_mumin) { feLogError("mumax must be > mumin"); return false; }
	return FEDamageCDF::Validate();
}

//-----------------------------------------------------------------------------
// Piecewise S-shaped quintic polynomial damage cumulative distribution function
double FEDamageCDFPQP::cdf(const double X)
{
    double cdf = 0;
    
    if (X <= m_mumin) cdf = 0;
    else if (X >= m_mumax) cdf = 1;
    else
    {
        double x = (X - m_mumin)/(m_mumax - m_mumin);
        cdf = pow(x,3)*(10 - 15*x + 6*x*x);
    }

    return cdf;
}

// Piecewise S-shaped quintic polynomial damage probability density function
double FEDamageCDFPQP::pdf(const double X)
{
    double pdf = 0;
    
    if (X <= m_mumin) pdf = 0;
    else if (X >= m_mumax) pdf = 0;
    else
    {
        double x = (X - m_mumin)/(m_mumax - m_mumin);
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
double FEDamageCDFGamma::cdf(const double X)
{
    double cdf = 0;
    
    // this CDF only admits positive values
#ifdef HAVE_GSL
    if (X > 0)
        cdf = gsl_sf_gamma_inc_P(m_alpha,X/m_mu);
#endif
    
    return cdf;
}

// Gamma damage probability density function
double FEDamageCDFGamma::pdf(const double X)
{
    double pdf = 0;
    
    // this CDF only admits positive values
#ifdef HAVE_GSL
    if (X > 0)
        pdf = pow(X/m_mu, m_alpha-1)*exp(-X/m_mu)/m_mu*gsl_sf_gammainv(m_alpha);
#endif
    
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
double FEDamageCDFUser::cdf(const double X)
{
    return m_cdf->value(X);
}

// Derivative of user-defined loadcurve damage probability density function
double FEDamageCDFUser::pdf(const double X)
{
    return m_cdf->derive(X);
}

