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
    ADD_PARAMETER(m_t[0], FE_PARAM_DOUBLE, "t1");
    ADD_PARAMETER(m_t[1], FE_PARAM_DOUBLE, "t2");
    ADD_PARAMETER(m_t[2], FE_PARAM_DOUBLE, "t3");
    ADD_PARAMETER(m_t[3], FE_PARAM_DOUBLE, "t4");
    ADD_PARAMETER(m_t[4], FE_PARAM_DOUBLE, "t5");
    ADD_PARAMETER(m_t[5], FE_PARAM_DOUBLE, "t6");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationExponential::FEBondRelaxationExponential(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_nt = 0;
	for (int i=0; i<MAX_TERMS; ++i) m_t[i] = 0;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEBondRelaxationExponential::Init()
{
    FEBondRelaxation::Init();
    
    m_nt = 0;
	for (int i=0; i<MAX_TERMS; ++i) {
        if (m_t[i] < 0) throw MaterialError("relaxation times must be > 0");
        if (m_t[i] > 0) ++m_nt;
    }
    if (m_nt == 0)
        throw MaterialError("at least one relaxation time must be > 0");
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationExponential::Relaxation(FEMaterialPoint& mp, const double t)
{
	// --- constant relaxation times ---
    double g = 0;
    
    for (int i=0; i<m_nt; ++i)
        if (m_t[i] > 0) {
            g += exp(-t/m_t[i]);
        }
	
	return g/m_nt;
}

///////////////////////////////////////////////////////////////////////////////
//
// FEBondRelaxationExpDistortion
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationExpDistortion, FEBondRelaxation)
ADD_PARAMETER(m_t[0], FE_PARAM_DOUBLE, "t1");
ADD_PARAMETER(m_t[1], FE_PARAM_DOUBLE, "t2");
ADD_PARAMETER(m_t[2], FE_PARAM_DOUBLE, "t3");
ADD_PARAMETER(m_t[3], FE_PARAM_DOUBLE, "t4");
ADD_PARAMETER(m_t[4], FE_PARAM_DOUBLE, "t5");
ADD_PARAMETER(m_t[5], FE_PARAM_DOUBLE, "t6");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationExpDistortion::FEBondRelaxationExpDistortion(FEModel* pfem) : FEBondRelaxation(pfem)
{
    m_nt = 0;
    for (int i=0; i<MAX_TERMS; ++i) m_t[i] = 0;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEBondRelaxationExpDistortion::Init()
{
    FEBondRelaxation::Init();
    
    m_nt = 0;
    for (int i=0; i<MAX_TERMS; ++i) {
        if (m_t[i] < 0) throw MaterialError("relaxation times must be > 0");
        if (m_t[i] > 0) ++m_nt;
    }
    if (m_nt == 0)
        throw MaterialError("at least one relaxation time must be > 0");
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationExpDistortion::Relaxation(FEMaterialPoint& mp, const double t)
{
    double eps = 1e-9;
    
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

    double g;
    
    // relax only if deformation is distortional
    if (K2 >= eps) {
        g = 0;
        for (int i=0; i<m_nt; ++i)
            if (m_t[i] > 0) {
                g += exp(-t/m_t[i]);
            }
        g /= m_nt;
    }
    else
        g = 1;
    
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
        g = (gsl_sf_expint_Ei(-t/m_tau1) - gsl_sf_expint_Ei(-t/m_tau2))/log(m_tau2/m_tau1);
    }
    else
        g = 1;
#endif
    
    return g;
}

///////////////////////////////////////////////////////////////////////////////
//
// FEBondRelaxationFungDistortional
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationFungDistortion, FEBondRelaxation)
ADD_PARAMETER(m_tau1, FE_PARAM_DOUBLE, "tau1");
ADD_PARAMETER(m_tau2, FE_PARAM_DOUBLE, "tau2");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationFungDistortion::FEBondRelaxationFungDistortion(FEModel* pfem) : FEBondRelaxation(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEBondRelaxationFungDistortion::Init()
{
    FEBondRelaxation::Init();
    
    if (m_tau1 <= 0) throw MaterialError("tau1 must be > 0");
    if (m_tau2 <= m_tau1) throw MaterialError("tau2 must be > tau1");
}

//-----------------------------------------------------------------------------
//! Relaxation function
double FEBondRelaxationFungDistortion::Relaxation(FEMaterialPoint& mp, const double t)
{
    double g = 0;
    
#ifdef HAVE_GSL
    double eps = 1e-9;
    
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
    
    // relax only if deformation is distortional
    if (K2 >= eps) {
        if (t > 0) {
            g = (gsl_sf_expint_Ei(-t/m_tau1) - gsl_sf_expint_Ei(-t/m_tau2))/log(m_tau2/m_tau1);
        }
        else
            g = 1;
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
// FEBondRelaxationFungDistortion
//
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEBondRelaxationParkDistortion, FEBondRelaxation)
ADD_PARAMETER(m_tau, FE_PARAM_DOUBLE, "tau");
ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEBondRelaxationParkDistortion::FEBondRelaxationParkDistortion(FEModel* pfem) : FEBondRelaxation(pfem)
{
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
    double eps = 1e-9;
    
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

    // relax only if deformation is distortional
    if (K2 >= eps)
        g = 1./(1+pow(t/m_tau,m_beta));
    else
        g = 1;
    
    return g;
}
