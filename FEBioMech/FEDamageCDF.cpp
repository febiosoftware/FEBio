//
//  FEDamageCDF.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/18/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEDamageCDF.h"
#include "FEDamageCriterion.h"
#include "FEDamageMaterialPoint.h"
#ifdef WIN32
#include "FECore\erf.h"
#endif

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEDamageCDFSimo, FEDamageCDF)
    ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "a");
    ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "b");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFSimo::FEDamageCDFSimo(FEModel* pfem) : FEDamageCDF(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEDamageCDFSimo::Init()
{
	if (m_alpha < 0) throw MaterialError("Invalid value of a: must be a non-negative number");
 	if (!INRANGE(m_beta, 0.0, 1.0)) throw MaterialError("Invalid value for b: must be in range [0,1]");
}

//-----------------------------------------------------------------------------
// Simo damage cumulative distribution function
// Simo, CMAME 60 (1987), 153-173
double FEDamageCDFSimo::Damage(FEMaterialPoint& mp)
{
    if (m_alpha == 0) {
        return 1;
    }
    
	// get the damage material point data
	FEDamageMaterialPoint& dp = *mp.ExtractData<FEDamageMaterialPoint>();
    
	// get the damage criterion (assuming dp.m_Etrial is up-to-date)
	double Es = max(dp.m_Etrial, dp.m_Emax);
    
    // evaluate the damage
    // this CDF only admits positive values
    if (Es >= 0) {
        double g = 1.0;
        if (fabs(Es) > 1e-12) g = m_beta + (1.0 - m_beta)*(1.0 - exp(-Es/m_alpha))/(Es/m_alpha);
        else g = 1.0 - 0.5*(1.0 - m_beta)/m_alpha*Es;
        dp.m_D = 1 - g;
    }
    
	return dp.m_D;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEDamageCDFLogNormal, FEDamageCDF)
    ADD_PARAMETER(m_mu, FE_PARAM_DOUBLE, "mu");
    ADD_PARAMETER(m_sigma, FE_PARAM_DOUBLE, "sigma");
    ADD_PARAMETER(m_Dmax, FE_PARAM_DOUBLE, "Dmax");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFLogNormal::FEDamageCDFLogNormal(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_mu = 1;
    m_sigma = 1;
    m_Dmax = 1;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEDamageCDFLogNormal::Init()
{
	if (m_mu <= 0) throw MaterialError("mu must be > 0");
	if (m_sigma <= 0) throw MaterialError("sigma must be > 0");
 	if (!INRANGE(m_Dmax, 0.0, 1.0)) throw MaterialError("Dmax must be in the range [0,1]");
}

//-----------------------------------------------------------------------------
//! Damage
double FEDamageCDFLogNormal::Damage(FEMaterialPoint& mp)
{
	// get the damage material point data
	FEDamageMaterialPoint& dp = *mp.ExtractData<FEDamageMaterialPoint>();
    
	// get the damage criterion (assuming dp.m_Etrial is up-to-date)
	double Es = max(dp.m_Etrial, dp.m_Emax);
    
    // evaluate the damage
    // this CDF only admits positive values
    if (Es > 0)
        dp.m_D = 0.5*erfc(-log(Es/m_mu)/m_sigma/sqrt(2.))*m_Dmax;
    
	return dp.m_D;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEDamageCDFWeibull, FEDamageCDF)
    ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
    ADD_PARAMETER(m_mu, FE_PARAM_DOUBLE, "mu");
    ADD_PARAMETER(m_Dmax, FE_PARAM_DOUBLE, "Dmax");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFWeibull::FEDamageCDFWeibull(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_alpha = m_mu = m_Dmax = 1;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEDamageCDFWeibull::Init()
{
	if (m_alpha <= 0) throw MaterialError("alpha must be >= 0");
	if (m_mu < 0) throw MaterialError("mu must be > 0");
 	if (!INRANGE(m_Dmax, 0.0, 1.0)) throw MaterialError("Dmax must be in the range [0,1]");
}

//-----------------------------------------------------------------------------
//! Damage
double FEDamageCDFWeibull::Damage(FEMaterialPoint& mp)
{
	// get the damage material point data
	FEDamageMaterialPoint& dp = *mp.ExtractData<FEDamageMaterialPoint>();
    
	// get the damage criterion (assuming dp.m_Etrial is up-to-date)
	double Es = max(dp.m_Etrial, dp.m_Emax);
    
    // evaluate the damage
    // this CDF only admits positive values
    if (Es > 0)
        dp.m_D = (1 - exp(-pow(Es/m_mu,m_alpha)))*m_Dmax;
    
	return dp.m_D;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEDamageCDFStep, FEDamageCDF)
    ADD_PARAMETER(m_mu, FE_PARAM_DOUBLE, "mu");
    ADD_PARAMETER(m_Dmax, FE_PARAM_DOUBLE, "Dmax");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFStep::FEDamageCDFStep(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_mu = m_Dmax = 1;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEDamageCDFStep::Init()
{
	if (m_mu < 0) throw MaterialError("mu must be >= 0");
 	if (!INRANGE(m_Dmax, 0.0, 1.0)) throw MaterialError("Dmax must be in the range [0,1]");
}

//-----------------------------------------------------------------------------
// Step cumulative distribution function (sudden fracture)
double FEDamageCDFStep::Damage(FEMaterialPoint& mp)
{
	// get the damage material point data
	FEDamageMaterialPoint& dp = *mp.ExtractData<FEDamageMaterialPoint>();
    
	// get the damage criterion (assuming dp.m_Etrial is up-to-date)
	double Es = max(dp.m_Etrial, dp.m_Emax);
    
    // evaluate the damage
    if (Es > m_mu) dp.m_D = m_Dmax;
    
	return dp.m_D;
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEDamageCDFPQP, FEDamageCDF)
    ADD_PARAMETER(m_mumin, FE_PARAM_DOUBLE, "mumin");
    ADD_PARAMETER(m_mumax, FE_PARAM_DOUBLE, "mumax");
    ADD_PARAMETER(m_Dmax, FE_PARAM_DOUBLE, "Dmax");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFPQP::FEDamageCDFPQP(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_mumin = 0;
    m_mumax = 1;
    m_Dmax = 1;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEDamageCDFPQP::Init()
{
	if (m_mumin < 0) throw MaterialError("mumin must be >= 0");
	if (m_mumax <= m_mumin) throw MaterialError("mumax must be > mumin");
 	if (!INRANGE(m_Dmax, 0.0, 1.0)) throw MaterialError("Dmax must be in the range [0,1]");
}

//-----------------------------------------------------------------------------
// Piecewise S-shaped quintic polynomial cumulative distribution function
double FEDamageCDFPQP::Damage(FEMaterialPoint& mp)
{
	// get the damage material point data
	FEDamageMaterialPoint& dp = *mp.ExtractData<FEDamageMaterialPoint>();
    
	// get the damage criterion (assuming dp.m_Etrial is up-to-date)
	double Es = max(dp.m_Etrial, dp.m_Emax);
    
    // evaluate the damage
	if (Es <= m_mumin) dp.m_D = 0;
	else if (Es >= m_mumax) dp.m_D = m_Dmax;
	else
	{
		double x = (Es - m_mumin)/(m_mumax - m_mumin);
        dp.m_D = pow(x,3)*(10 - 15*x + 6*x*x)*m_Dmax;
	}
    
	return dp.m_D;
}
