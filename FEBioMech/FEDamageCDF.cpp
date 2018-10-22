//
//  FEDamageCDF.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/18/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEDamageCDF.h"
#include "FEDamageCriterion.h"
#include "FEDamageMaterialPoint.h"
#define _USE_MATH_DEFINES
#include <math.h>

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
    ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(1.0), "alpha");
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
	if (m_mumax <= m_mumin) return fecore_error("mumax must be > mumin");
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
