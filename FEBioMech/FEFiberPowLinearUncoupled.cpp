//
//  FEFiberPowLinearUncoupled.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/6/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FEFiberPowLinearUncoupled.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberPowLinearUncoupled, FEUncoupledMaterial)
	ADD_PARAMETER2(m_E    , FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "E"    );
	ADD_PARAMETER2(m_lam0 , FE_PARAM_DOUBLE, FE_RANGE_GREATER(1.0), "lam0" );
	ADD_PARAMETER2(m_beta , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta" );
	ADD_PARAMETER(m_thd  , FE_PARAM_DOUBLE, "theta");
	ADD_PARAMETER(m_phd  , FE_PARAM_DOUBLE, "phi"  );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEFiberPowLinearUncoupled.h
//-----------------------------------------------------------------------------

void FEFiberPowLinearUncoupled::Init()
{
	FEElasticMaterial::Init();

    // convert angles from degrees to radians
    double pi = 4*atan(1.0);
    double the = m_thd*pi/180.;
    double phi = m_phd*pi/180.;
    // fiber direction in local coordinate system (reference configuration)
    m_n0.x = cos(the)*sin(phi);
    m_n0.y = sin(the)*sin(phi);
    m_n0.z = cos(phi);
    
    // initialize material constants
    m_I0 = m_lam0*m_lam0;
    m_ksi = m_E/(4*m_I0*m_I0)/(m_beta-1)*pow(m_I0-1, 2-m_beta);
    m_b = m_ksi*pow(m_I0-1, m_beta-1) + m_E/4/m_I0;
}

//-----------------------------------------------------------------------------
mat3ds FEFiberPowLinearUncoupled::DevStress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    double J = pt.m_J;
    mat3d F = pt.m_F*pow(J, -1./3.);
    
    // loop over all integration points
    vec3d n0, nt;
    double In, sn;
    const double eps = 0;
    mat3ds C = pt.DevRightCauchyGreen();
    mat3ds s;
    
    // evaluate fiber direction in global coordinate system
    n0 = pt.m_Q*m_n0;
    
    // Calculate In
    In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 > eps)
    {
        // get the global spatial fiber direction in current configuration
        nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate the fiber stress magnitude
        sn = (In < m_I0) ?
        2*In*m_ksi*pow(In-1, m_beta-1) :
        2*m_b*In - m_E/2;
        
        // calculate the fiber stress
        s = N*(sn/J);
    }
    else
    {
        s.zero();
    }
    
    return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEFiberPowLinearUncoupled::DevTangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    double J = pt.m_J;
    mat3d F = pt.m_F*pow(J, -1./3.);
    
    // loop over all integration points
    vec3d n0, nt;
    double In, sn, cn;
    const double eps = 0;
    mat3ds C = pt.DevRightCauchyGreen();
    mat3ds s;
    tens4ds c;
    
    // evaluate fiber direction in global coordinate system
    n0 = pt.m_Q*m_n0;
    
    // Calculate In
    In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 > eps)
    {
        // get the global spatial fiber direction in current configuration
        nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // calculate the fiber stress magnitude
        sn = (In < m_I0) ?
        2*In*m_ksi*pow(In-1, m_beta-1) :
        2*m_b*In - m_E/2;
        
        // calculate the fiber stress
        s = N*(sn/J);
        
        // calculate modulus
        cn = (In < m_I0) ?
        4*In*In*m_ksi*(m_beta-1)*pow(In-1, m_beta-2) :
        m_E;
        
        // calculate the fiber tangent
        c = NxN*(cn/J);

        // This is the final value of the elasticity tensor
        mat3dd I(1);
        tens4ds IxI = dyad1s(I);
        tens4ds I4  = dyad4s(I);
        c += ((I4+IxI/3.0)*s.tr() - dyad1s(I,s))*(2./3.)
        - (ddots(IxI, c)-IxI*(c.tr()/3.))/3.;
    }
    else
    {
        c.zero();
    }
    
    return c;
}

//-----------------------------------------------------------------------------
double FEFiberPowLinearUncoupled::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // loop over all integration points
    vec3d n0;
    double In;
    const double eps = 0;
    mat3ds C = pt.DevRightCauchyGreen();
    
    // evaluate fiber direction in global coordinate system
    n0 = pt.m_Q*m_n0;
    
    // Calculate In = n0*C*n0
    In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 > eps)
    {
        // calculate strain energy density
        sed = (In < m_I0) ?
        m_ksi/m_beta*pow(In-1, m_beta) :
        m_b*(In-m_I0) - m_E/4*log(In/m_I0) + m_ksi/m_beta*pow(m_I0-1, m_beta);
    }
    
    return sed;
}
