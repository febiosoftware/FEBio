//
//  FEFiberIntegrationGeodesicUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//
/*
#include "FEFiberIntegrationGeodesicUC.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FEFiberIntegrationGeodesicUC
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberIntegrationGeodesicUC, FEFiberIntegrationSchemeUC)
    ADD_PARAMETER(m_nres, FE_PARAM_INT, "resolution");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEFiberIntegrationGeodesicUC::Serialize(DumpStream& ar)
{
    FEFiberIntegrationSchemeUC::Serialize(ar);
    if (ar.IsSaving() == false)
    {
        InitIntegrationRule();
    }
}

//-----------------------------------------------------------------------------
bool FEFiberIntegrationGeodesicUC::Init()
{
	if ((m_nres != 0) && (m_nres != 1)) return MaterialError("resolution must be 0 (low) or 1 (high).");
    
    // initialize integration rule data
    InitIntegrationRule();
    
    // also initialize the parent class
    return FEFiberIntegrationSchemeUC::Init();
}

//-----------------------------------------------------------------------------
void FEFiberIntegrationGeodesicUC::InitIntegrationRule()
{
    // select the integration rule
    m_nint = (m_nres == 0? NSTL  : NSTH  );
    const double* phi = (m_nres == 0? PHIL  : PHIH  );
    const double* the = (m_nres == 0? THETAL: THETAH);
    const double* w   = (m_nres == 0? AREAL : AREAH );
    
    for (int n=0; n<m_nint; ++n)
    {
        m_cth[n] = cos(the[n]);
        m_sth[n] = sin(the[n]);
        m_cph[n] = cos(phi[n]);
        m_sph[n] = sin(phi[n]);
        m_w[n] = w[n];
    }
}

//-----------------------------------------------------------------------------
mat3ds FEFiberIntegrationGeodesicUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();
    
	// loop over all integration points
	double R;
	vec3d n0e, n0a;
	mat3ds s;
	s.zero();
    
	for (int n=0; n<m_nint; ++n)
	{
		// set the global fiber direction in reference configuration
		n0e.x = m_cth[n]*m_sph[n];
		n0e.y = m_sth[n]*m_sph[n];
		n0e.z = m_cph[n];
        m_pFmat->SetFiberDirection(mp, n0e);
        
        // get the local material fiber direction in reference configuration
        n0a = QT*n0e;
        // evaluate the fiber density
        R = m_pFDD->FiberDensity(n0a);
        
        // evaluate this fiber's contribution to the stress
		s += m_pFmat->DevStress(mp)*(R*m_w[n]);
	}
    // we don't need to evaluate the deviatoric part since we already summed up deviatoric fiber stresses
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberIntegrationGeodesicUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();
    
	// loop over all integration points
	double R;
	vec3d n0e, n0a;
	tens4ds c;
	c.zero();
	
	for (int n=0; n<m_nint; ++n)
	{
		// set the global fiber direction in reference configuration
		n0e.x = m_cth[n]*m_sph[n];
		n0e.y = m_sth[n]*m_sph[n];
		n0e.z = m_cph[n];
        m_pFmat->SetFiberDirection(mp, n0e);
        
        // get the local material fiber direction in reference configuration
        n0a = QT*n0e;
        // evaluate the fiber density
        R = m_pFDD->FiberDensity(n0a);
        
        // evaluate this fiber's contribution to the tangent
		c += m_pFmat->DevTangent(mp)*(R*m_w[n]);
	}
	
    // we don't need to evaluate the deviatoric part since we already summed up deviatoric fiber tangents
	return c;
}

//-----------------------------------------------------------------------------
double FEFiberIntegrationGeodesicUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();
    
	// loop over all integration points
	double R;
	vec3d n0e, n0a;
	double sed = 0.0;
    
	for (int n=0; n<m_nint; ++n)
	{
		// set the global fiber direction in reference configuration
		n0e.x = m_cth[n]*m_sph[n];
		n0e.y = m_sth[n]*m_sph[n];
		n0e.z = m_cph[n];
        m_pFmat->SetFiberDirection(mp, n0e);
        
        // get the local material fiber direction in reference configuration
        n0a = QT*n0e;
        // evaluate the fiber density
        R = m_pFDD->FiberDensity(n0a);
        
        // evaluate this fiber's contribution to the stress
		sed += m_pFmat->DevStrainEnergyDensity(mp)*(R*m_w[n]);
	}
    
	return sed;
}

//-----------------------------------------------------------------------------
double FEFiberIntegrationGeodesicUC::IntegratedFiberDensity()
{
    // initialize integrated fiber density distribution
    double IFD = 1;
    
	// loop over all integration points
	double R;
	vec3d n0a;
    double C = 0;
    
	for (int n=0; n<m_nint; ++n)
	{
		// set the global fiber direction in reference configuration
		n0a.x = m_cth[n]*m_sph[n];
		n0a.y = m_sth[n]*m_sph[n];
		n0a.z = m_cph[n];
        
        // evaluate the fiber density
        R = m_pFDD->FiberDensity(n0a);
        
        // integrate the fiber density
		C += R*m_w[n];
	}
	IFD = C;
    return IFD;
}
*/
