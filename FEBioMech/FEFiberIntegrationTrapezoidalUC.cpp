//
//  FEFiberIntegrationTrapezoidalUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEFiberIntegrationTrapezoidalUC.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FEFiberIntegrationTrapezoidalUC
//-----------------------------------------------------------------------------

// register the material with the framework

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberIntegrationTrapezoidalUC, FEFiberIntegrationSchemeUC)
    ADD_PARAMETER(m_nth, FE_PARAM_INT, "nth");
END_PARAMETER_LIST();

void FEFiberIntegrationTrapezoidalUC::Init()
{
	if (m_nth < 1) throw MaterialError("nth must be strictly greater than zero.");
    
    // also initialize the parent class
    FEFiberIntegrationSchemeUC::Init();
}

//-----------------------------------------------------------------------------
mat3ds FEFiberIntegrationTrapezoidalUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
    // initialize stress tensor
	mat3ds s;
	s.zero();
    
    double theta;
    double pi = 4*atan(1.0);
    double dth = pi/m_nth;  // integrate from 0 to pi
    
    // get the element's local coordinate system
	mat3d Q = pt.m_Q;
    vec3d a0(Q(0,0),Q(1,0),Q(2,0)); // local x-direction unit vector
    vec3d a1(Q(0,1),Q(1,1),Q(2,1)); // local y-direction unit vector
    
    vec3d n0e, n0a;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        
        // set fiber direction in x-y plane of local coordinate system
        n0a = a0*cos(theta) + a1*sin(theta);
        // evaluate local fiber distribution
        double R = m_pFDD->FiberDensity(n0a);
        
        // rotate to global configuration to set fiber direction
        n0e = Q*n0a;
        m_pFmat->SetFiberDirection(mp, n0e);
        
        // calculate the stress
        s += m_pFmat->DevStress(mp)*(R*dth);
    }
    
    // Multiply by 2 since fibers along theta+pi have same stress as along theta
    // we don't need to evaluate the deviatoric part since we already summed up deviatoric fiber stresses
	return s*2;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberIntegrationTrapezoidalUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
    // initialize stress tensor
	tens4ds c;
	c.zero();
    
    double theta;
    double pi = 4*atan(1.0);
    double dth = pi/m_nth;  // integrate from 0 to pi
    
    // get the element's local coordinate system
	mat3d Q = pt.m_Q;
    vec3d a0(Q(0,0),Q(1,0),Q(2,0)); // local x-direction unit vector
    vec3d a1(Q(0,1),Q(1,1),Q(2,1)); // local y-direction unit vector
    
    vec3d n0e, n0a;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        
        // set fiber direction in x-y plane of local coordinate system
        n0a = a0*cos(theta) + a1*sin(theta);
        // evaluate local fiber distribution
        double R = m_pFDD->FiberDensity(n0a);
        
        // rotate to global configuration to set fiber direction
        n0e = Q*n0a;
        m_pFmat->SetFiberDirection(mp, n0e);
        
        // calculate the stress
        c += m_pFmat->DevTangent(mp)*(R*dth);
    }
    
    // Multiply by 2 since fibers along theta+pi have same stress as along theta
    // we don't need to evaluate the deviatoric part since we already summed up deviatoric fiber tangents
    return c*2;
}

//-----------------------------------------------------------------------------
double FEFiberIntegrationTrapezoidalUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
    // initialize strain energy density
	double sed = 0.0;
    
    double theta;
    double pi = 4*atan(1.0);
    double dth = pi/m_nth;  // integrate from 0 to pi
    
    // get the element's local coordinate system
	mat3d Q = pt.m_Q;
    vec3d a0(Q(0,0),Q(1,0),Q(2,0)); // local x-direction unit vector
    vec3d a1(Q(0,1),Q(1,1),Q(2,1)); // local y-direction unit vector
    
    vec3d n0e, n0a;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        
        // set fiber direction in x-y plane of local coordinate system
        n0a = a0*cos(theta) + a1*sin(theta);
        // evaluate local fiber distribution
        double R = m_pFDD->FiberDensity(n0a);
        
        // rotate to global configuration to set fiber direction
        n0e = Q*n0a;
        m_pFmat->SetFiberDirection(mp, n0e);
        
        // calculate the stress
        sed += m_pFmat->DevStrainEnergyDensity(mp)*(R*dth);
    }
    
    // Multiply by 2 since fibers along theta+pi have same sed as along theta
	return sed*2;
}

//-----------------------------------------------------------------------------
void FEFiberIntegrationTrapezoidalUC::IntegratedFiberDensity(double& IFD)
{
    // initialize integrated fiber density distribution
    IFD = 1;
    
	double C = 0;
    
    double theta;
    double pi = 4*atan(1.0);
    double dth = pi/m_nth;  // integrate from 0 to pi
    
    vec3d a0(1,0,0); // local x-direction unit vector
    vec3d a1(0,1,0); // local y-direction unit vector
    
    vec3d n0a;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        
        // set fiber direction in x-y plane of local coordinate system
        n0a = a0*cos(theta) + a1*sin(theta);
        // evaluate local fiber distribution
        double R = m_pFDD->FiberDensity(n0a);
        
        // integrate the fiber distribution
        C += R*dth;
    }
    
    // Multiply by 2 to take advantange of symmetry
	IFD = C*2;
}
