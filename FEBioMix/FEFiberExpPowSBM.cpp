//
//  FEFiberExpPowSBM.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 5/24/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FEFiberExpPowSBM.h"
#include "FEMultiphasic.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberExpPowSBM, FEElasticMaterial)
ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
ADD_PARAMETER(m_beta , FE_PARAM_DOUBLE, "beta" );
ADD_PARAMETER(m_ksi0 , FE_PARAM_DOUBLE, "ksi0" );
ADD_PARAMETER(m_rho0 , FE_PARAM_DOUBLE, "rho0" );
ADD_PARAMETER(m_g    , FE_PARAM_DOUBLE, "gamma");
ADD_PARAMETER(m_sbm  , FE_PARAM_INT   , "sbm"  );
ADD_PARAMETER(m_thd  , FE_PARAM_DOUBLE, "theta");
ADD_PARAMETER(m_phd  , FE_PARAM_DOUBLE, "phi"  );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEFiberExpPow
//-----------------------------------------------------------------------------

void FEFiberExpPowSBM::Init()
{
    if (m_ksi0 < 0) throw MaterialError("ksi0 must be positive.");
    if (m_beta < 2) throw MaterialError("beta must be >= 2.");
    if (m_alpha < 0) throw MaterialError("alpha must be >= 0.");
    if (m_rho0 < 0) throw MaterialError("rho0 must be positive.");
    if (m_g < 0) throw MaterialError("gamma must be positive.");
    
    // get the parent material which must be a multiphasic material
    FEMultiphasic* pMP = dynamic_cast<FEMultiphasic*> (GetAncestor());
    if (pMP == 0) throw MaterialError("Parent material must be multiphasic");
    
    // extract the local id of the SBM whose density controls Young's modulus from the global id
    m_lsbm = pMP->FindLocalSBMID(m_sbm);
    if (m_lsbm == -1) throw MaterialError("Invalid value for sbm");
    
    
    // convert angles from degrees to radians
    double pi = 4*atan(1.0);
    double the = m_thd*pi/180.;
    double phi = m_phd*pi/180.;
    // fiber direction in local coordinate system (reference configuration)
    m_n0.x = cos(the)*sin(phi);
    m_n0.y = sin(the)*sin(phi);
    m_n0.z = cos(phi);
}

//-----------------------------------------------------------------------------
mat3ds FEFiberExpPowSBM::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    // initialize material constants
    double rhor = spt.m_sbmr[m_lsbm];
    double ksi = FiberModulus(rhor);
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    vec3d n0, nt;
    double In_1, Wl;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    mat3ds s;
    
    // evaluate fiber direction in global coordinate system
    n0 = pt.m_Q*m_n0;
    
    // Calculate In = n0*C*n0
    In_1 = n0*(C*n0) - 1.0;
    
    // only take fibers in tension into consideration
    if (In_1 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        nt = F*n0;
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate strain energy derivative
        Wl = ksi*pow(In_1, m_beta-1.0)*exp(m_alpha*pow(In_1, m_beta));
        
        // calculate the fiber stress
        s = N*(2.0*Wl/J);
    }
    else
    {
        s.zero();
    }
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberExpPowSBM::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    // initialize material constants
    double rhor = spt.m_sbmr[m_lsbm];
    double ksi = FiberModulus(rhor);
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    vec3d n0, nt;
    double In_1, Wll;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    tens4ds c;
    
    // evaluate fiber direction in global coordinate system
    n0 = pt.m_Q*m_n0;
    
    // Calculate In = n0*C*n0
    In_1 = n0*(C*n0) - 1.0;
    
    // only take fibers in tension into consideration
    if (In_1 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        nt = F*n0;
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // calculate strain energy 2nd derivative
        double tmp = m_alpha*pow(In_1, m_beta);
        Wll = ksi*pow(In_1, m_beta-2.0)*((tmp+1)*m_beta-1.0)*exp(tmp);
        
        // calculate the fiber tangent
        c = NxN*(4.0*Wll/J);
    }
    else
    {
        c.zero();
    }
    
    return c;
}

//-----------------------------------------------------------------------------
double FEFiberExpPowSBM::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    // initialize material constants
    double rhor = spt.m_sbmr[m_lsbm];
    double ksi = FiberModulus(rhor);
    
    // loop over all integration points
    vec3d n0;
    double In_1;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    
    // evaluate fiber direction in global coordinate system
    n0 = pt.m_Q*m_n0;
    
    // Calculate In = n0*C*n0
    In_1 = n0*(C*n0) - 1.0;
    
    // only take fibers in tension into consideration
    if (In_1 >= eps)
    {
        // calculate strain energy derivative
        if (m_alpha > 0) {
            sed = ksi/(m_alpha*m_beta)*(exp(m_alpha*pow(In_1, m_beta))-1);
        }
        else
            sed = ksi/m_beta*pow(In_1, m_beta);
    }
    
    return sed;
}
