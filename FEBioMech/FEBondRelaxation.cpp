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

//-----------------------------------------------------------------------------
// Material parameters for FEBondRelaxation
void FEBondRelaxation::Init()
{
	FEMaterial::Init();
}

///////////////////////////////////////////////////////////////////////////////
//
// FEBondRelaxationExponential
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationExponential, FEBondRelaxation)
    ADD_PARAMETER(m_tau, FE_PARAM_DOUBLE, "tau");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationExponential::FEBondRelaxationExponential(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_tau = 0;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEBondRelaxationExponential::Init()
{
    FEBondRelaxation::Init();
    
    if (m_tau <= 0) throw MaterialError("tau must be > 0");
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationExponential::Relaxation(FEMaterialPoint& mp, const double t)
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
    ADD_PARAMETER(m_tau, FE_PARAM_DOUBLE, "tau");
    ADD_PARAMETER(m_tau1, FE_PARAM_DOUBLE, "tau1");
    ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationExpDistortion::FEBondRelaxationExpDistortion(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_tau = 0;
    m_tau1 = 0;
    m_alpha = 1;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEBondRelaxationExpDistortion::Init()
{
    FEBondRelaxation::Init();
    
    if (m_tau <= 0) throw MaterialError("tau must be > 0");
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationExpDistortion::Relaxation(FEMaterialPoint& mp, const double t)
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
    double tau = m_tau + m_tau1*K2a;
    
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
    ADD_PARAMETER(m_tau1, FE_PARAM_DOUBLE, "tau1");
    ADD_PARAMETER(m_tau2, FE_PARAM_DOUBLE, "tau2");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationFung::FEBondRelaxationFung(FEModel* pfem) : FEBondRelaxation(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEBondRelaxationFung::Init()
{
    FEBondRelaxation::Init();
    
    if (m_tau1 <= 0) throw MaterialError("tau1 must be > 0");
    if (m_tau2 <= m_tau1) throw MaterialError("tau2 must be > tau1");
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationFung::Relaxation(FEMaterialPoint& mp, const double t)
{
    double g = 0;
    
#ifdef HAVE_GSL
    if (t > 0) {
        g = (m_tau2*exp(-t/m_tau2) - m_tau1*exp(-t/m_tau1)
        + t*(gsl_sf_expint_Ei(-t/m_tau1) - gsl_sf_expint_Ei(-t/m_tau2)))
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
    ADD_PARAMETER(m_tau, FE_PARAM_DOUBLE, "tau");
    ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationPark::FEBondRelaxationPark(FEModel* pfem) : FEBondRelaxation(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEBondRelaxationPark::Init()
{
    FEBondRelaxation::Init();
    
    if (m_tau <= 0) throw MaterialError("tau must be > 0");
    if (m_beta <= 0) throw MaterialError("beta must be > 0");
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationPark::Relaxation(FEMaterialPoint& mp, const double t)
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
    ADD_PARAMETER(m_tau, FE_PARAM_DOUBLE, "tau");
    ADD_PARAMETER(m_tau1, FE_PARAM_DOUBLE, "tau1");
    ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
    ADD_PARAMETER(m_beta1, FE_PARAM_DOUBLE, "beta1");
    ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationParkDistortion::FEBondRelaxationParkDistortion(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_tau1 = 0;
    m_beta1 = 0;
    m_alpha = 1;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEBondRelaxationParkDistortion::Init()
{
    FEBondRelaxation::Init();
    
    if (m_tau <= 0) throw MaterialError("tau must be > 0");
    if (m_beta <= 0) throw MaterialError("beta must be > 0");
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationParkDistortion::Relaxation(FEMaterialPoint& mp, const double t)
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
    double tau = m_tau + m_tau1*K2a;
    double beta = m_beta + m_beta1*K2a;
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
    ADD_PARAMETER(m_tau, FE_PARAM_DOUBLE, "tau");
    ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationPower::FEBondRelaxationPower(FEModel* pfem) : FEBondRelaxation(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEBondRelaxationPower::Init()
{
    FEBondRelaxation::Init();
    
    if (m_tau <= 0) throw MaterialError("tau must be > 0");
    if (m_beta <= 0) throw MaterialError("beta must be > 0");
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationPower::Relaxation(FEMaterialPoint& mp, const double t)
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
    ADD_PARAMETER(m_tau, FE_PARAM_DOUBLE, "tau");
    ADD_PARAMETER(m_tau1, FE_PARAM_DOUBLE, "tau1");
    ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
    ADD_PARAMETER(m_beta1, FE_PARAM_DOUBLE, "beta1");
    ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationPowerDistortion::FEBondRelaxationPowerDistortion(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_tau1 = 0;
    m_alpha = 1;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEBondRelaxationPowerDistortion::Init()
{
    FEBondRelaxation::Init();
    
    if (m_tau <= 0) throw MaterialError("tau must be > 0");
    if (m_beta <= 0) throw MaterialError("beta must be > 0");
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationPowerDistortion::Relaxation(FEMaterialPoint& mp, const double t)
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
    double tau = m_tau + m_tau1*K2a;
    double beta = m_beta + m_beta1*K2a;

    g = pow(1+t/tau,-beta);
    
    return g;
}
