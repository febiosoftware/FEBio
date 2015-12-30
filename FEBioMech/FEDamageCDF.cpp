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
    ADD_PARAMETER2(m_alpha, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "a");
    ADD_PARAMETER2(m_beta, FE_PARAM_DOUBLE, FE_RANGE_CLOSED(0.0, 1.0), "b");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFSimo::FEDamageCDFSimo(FEModel* pfem) : FEDamageCDF(pfem)
{
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
    ADD_PARAMETER2(m_mu   , FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "mu");
    ADD_PARAMETER2(m_sigma, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "sigma");
    ADD_PARAMETER2(m_Dmax , FE_PARAM_DOUBLE, FE_RANGE_CLOSED(0.0, 1.0), "Dmax");
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
    ADD_PARAMETER2(m_alpha, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "alpha");
    ADD_PARAMETER2(m_mu   , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
    ADD_PARAMETER2(m_Dmax , FE_PARAM_DOUBLE, FE_RANGE_CLOSED(0.0, 1.0), "Dmax");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFWeibull::FEDamageCDFWeibull(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_alpha = m_mu = m_Dmax = 1;
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
    ADD_PARAMETER2(m_mu  , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
    ADD_PARAMETER2(m_Dmax, FE_PARAM_DOUBLE, FE_RANGE_CLOSED(0.0, 1.0), "Dmax");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEDamageCDFStep::FEDamageCDFStep(FEModel* pfem) : FEDamageCDF(pfem)
{
    m_mu = m_Dmax = 1;
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
    ADD_PARAMETER2(m_mumin, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mumin");
    ADD_PARAMETER(m_mumax , FE_PARAM_DOUBLE, "mumax");
    ADD_PARAMETER2(m_Dmax , FE_PARAM_DOUBLE, FE_RANGE_CLOSED(0.0, 1.0), "Dmax");
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
bool FEDamageCDFPQP::Validate()
{
	if (m_mumax <= m_mumin) return MaterialError("mumax must be > mumin");
	return FEDamageCDF::Validate();
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
