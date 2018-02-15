//
//  FEBondRelaxation.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/25/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEBondRelaxation.h"
#include "FEElasticMaterial.h"
#ifdef HAVE_GSL
    #include "gsl/gsl_sf_expint.h"
#endif

///////////////////////////////////////////////////////////////////////////////
//
// FEBondRelaxationExponential
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationExponential, FEBondRelaxation)
    ADD_PARAMETER2(m_tau, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "tau");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationExponential::FEBondRelaxationExponential(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_tau = 0;
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationExponential::Relaxation(FEMaterialPoint& mp, const double t, const mat3ds D)
{
	// --- constant relaxation times ---
    double g = exp(-t/m_tau);
	
	return g;
}

///////////////////////////////////////////////////////////////////////////////
//
// FEBondRelaxationExpDistortion
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationExpDistortion, FEBondRelaxation)
    ADD_PARAMETER2(m_tau0 , FE_PARAM_DOUBLE, FE_RANGE_GREATER         (0.0), "tau0" );
    ADD_PARAMETER2(m_tau1 , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "tau1" );
    ADD_PARAMETER2(m_alpha, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationExpDistortion::FEBondRelaxationExpDistortion(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_tau0 = 0;
    m_tau1 = 0;
    m_alpha = 1;
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationExpDistortion::Relaxation(FEMaterialPoint& mp, const double t, const mat3ds D)
{
    // get the elastic material point data
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // evaluate spatial Hencky (logarithmic) strain
    mat3ds B = pt.LeftCauchyGreen();
    double d[3];
    vec3d v[3];
    B.eigen2(d,v);
    mat3ds h = (dyad(v[0])*log(d[0]) + dyad(v[1])*log(d[1]) + dyad(v[2])*log(d[2]))/2;
    
    // evaluate distortion magnitude (always positive)
    double K2 = (h.dev()).norm();

    double K2a = pow(K2,m_alpha);
    double tau = m_tau0 + m_tau1*K2a;
    
    double g = exp(-t/tau);
    
    return g;
}

///////////////////////////////////////////////////////////////////////////////
//
// FEBondRelaxationFung
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationFung, FEBondRelaxation)
    ADD_PARAMETER2(m_tau1, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "tau1");
    ADD_PARAMETER2(m_tau2, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "tau2");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationFung::FEBondRelaxationFung(FEModel* pfem) : FEBondRelaxation(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEBondRelaxationFung::Validate()
{
    if (FEBondRelaxation::Validate() == false) return false;
    if (m_tau2 <= m_tau1) return MaterialError("tau2 must be > tau1");
	return true;
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationFung::Relaxation(FEMaterialPoint& mp, const double t, const mat3ds D)
{
    double g = 0;
    
#ifdef HAVE_GSL
    if (t > 0) {
        g = (m_tau2*exp(-t/m_tau2) - m_tau1*exp(-t/m_tau1)
        + t*(gsl_sf_expint_Ei(-t/m_tau2) - gsl_sf_expint_Ei(-t/m_tau1)))
        /(m_tau2 - m_tau1);
    }
    else
        g = 1;
#endif
    
    return g;
}

///////////////////////////////////////////////////////////////////////////////
//
// FEBondRelaxationPark
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationPark, FEBondRelaxation)
    ADD_PARAMETER2(m_tau , FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "tau");
    ADD_PARAMETER2(m_beta, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "beta");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationPark::FEBondRelaxationPark(FEModel* pfem) : FEBondRelaxation(pfem)
{
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationPark::Relaxation(FEMaterialPoint& mp, const double t, const mat3ds D)
{
    double g = 1./(1+pow(t/m_tau,m_beta));
    
    return g;
}

///////////////////////////////////////////////////////////////////////////////
//
// FEBondRelaxationParkDistortion
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationParkDistortion, FEBondRelaxation)
    ADD_PARAMETER2(m_tau0 , FE_PARAM_DOUBLE, FE_RANGE_GREATER         (0.0), "tau0" );
    ADD_PARAMETER2(m_tau1 , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "tau1" );
    ADD_PARAMETER2(m_beta0, FE_PARAM_DOUBLE, FE_RANGE_GREATER         (0.0), "beta0");
    ADD_PARAMETER2(m_beta1, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "beta1");
    ADD_PARAMETER2(m_alpha, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationParkDistortion::FEBondRelaxationParkDistortion(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_tau0 = 0;
    m_beta0 = 0;
    m_tau1 = 0;
    m_beta1 = 0;
    m_alpha = 1;
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationParkDistortion::Relaxation(FEMaterialPoint& mp, const double t, const mat3ds D)
{
    double g;
    
    // get the elastic material point data
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // evaluate spatial Hencky (logarithmic) strain
    mat3ds B = pt.LeftCauchyGreen();
    double d[3];
    vec3d v[3];
    B.eigen2(d,v);
    mat3ds h = (dyad(v[0])*log(d[0]) + dyad(v[1])*log(d[1]) + dyad(v[2])*log(d[2]))/2;
    
    // evaluate distortion magnitude (always positive)
    double K2 = (h.dev()).norm();

    double K2a = pow(K2,m_alpha);
    double tau = m_tau0 + m_tau1*K2a;
    double beta = m_beta0 + m_beta1*K2a;
    g = 1./(1+pow(t/tau,beta));
    
    return g;
}

///////////////////////////////////////////////////////////////////////////////
//
// FEBondRelaxationPower
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationPower, FEBondRelaxation)
    ADD_PARAMETER2(m_tau , FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "tau");
    ADD_PARAMETER2(m_beta, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "beta");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationPower::FEBondRelaxationPower(FEModel* pfem) : FEBondRelaxation(pfem)
{
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationPower::Relaxation(FEMaterialPoint& mp, const double t, const mat3ds D)
{
    double g = pow(1+t/m_tau,-m_beta);
    
    return g;
}

///////////////////////////////////////////////////////////////////////////////
//
// FEBondRelaxationPowerDistortion
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationPowerDistortion, FEBondRelaxation)
    ADD_PARAMETER2(m_tau0 , FE_PARAM_DOUBLE, FE_RANGE_GREATER         (0.0), "tau0");
    ADD_PARAMETER2(m_tau1 , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "tau1");
    ADD_PARAMETER2(m_beta0, FE_PARAM_DOUBLE, FE_RANGE_GREATER         (0.0), "beta0");
    ADD_PARAMETER2(m_beta1, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "beta1");
    ADD_PARAMETER2(m_alpha, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationPowerDistortion::FEBondRelaxationPowerDistortion(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_tau1 = 0;
    m_alpha = 1;
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationPowerDistortion::Relaxation(FEMaterialPoint& mp, const double t, const mat3ds D)
{
    double g;
    
    // get the elastic material point data
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // evaluate spatial Hencky (logarithmic) strain
    mat3ds B = pt.LeftCauchyGreen();
    double d[3];
    vec3d v[3];
    B.eigen2(d,v);
    mat3ds h = (dyad(v[0])*log(d[0]) + dyad(v[1])*log(d[1]) + dyad(v[2])*log(d[2]))/2;
    
    // evaluate distortion magnitude (always positive)
    double K2 = (h.dev()).norm();
    
    double K2a = pow(K2,m_alpha);
    double tau = m_tau0 + m_tau1*K2a;
    double beta = m_beta0 + m_beta1*K2a;

    g = pow(1+t/tau,-beta);
    
    return g;
}

///////////////////////////////////////////////////////////////////////////////
//
// FEBondRelaxationCarreau
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationCarreau, FEBondRelaxation)
    ADD_PARAMETER2(m_tau0 , FE_PARAM_DOUBLE, FE_RANGE_GREATER         (0.0), "tau0");
    ADD_PARAMETER2(m_lam  , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "lambda");
    ADD_PARAMETER2(m_n    , FE_PARAM_DOUBLE, FE_RANGE_GREATER         (0.0), "n");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationCarreau::FEBondRelaxationCarreau(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_tau0 = 0;
    m_lam = 0;
    m_n = 1;
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationCarreau::Relaxation(FEMaterialPoint& mp, const double t, const mat3ds D)
{
    double g;
    
    // evaluate the engineering shear rate
    double gdot = sqrt(2.)*D.norm();
    
    // evaluate the relaxation time
    double tau = m_tau0*pow(1+pow(m_lam*gdot,2),(m_n-1)/2.);
    
    g = exp(-t/tau);
    
    return g;
}
