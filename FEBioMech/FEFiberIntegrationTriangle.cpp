//
//  FEFiberIntegrationTriangle.cpp
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 11/19/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#include "stdafx.h"
#include "FEFiberIntegrationTriangle.h"
#include "triangle_sphere.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FEFiberIntegrationTriangle
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberIntegrationTriangle, FEFiberIntegrationScheme)
    ADD_PARAMETER(m_nres, FE_PARAM_INT, "resolution");
END_PARAMETER_LIST();

void FEFiberIntegrationTriangle::Init()
{
	m_bfirst = true;
	
	if (m_bfirst)
	{
		switch (m_nres) {
            
            
            case 1:
                m_nint=NST1;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA1[n]);
                    m_sth[n] = sin(THETA1[n]);
                    m_cph[n] = cos(PHI1[n]);
                    m_sph[n] = sin(PHI1[n]);
                    m_w[n] = AREA1[n];
                }
                break;
            case 2:
                m_nint=NST2;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA2[n]);
                    m_sth[n] = sin(THETA2[n]);
                    m_cph[n] = cos(PHI2[n]);
                    m_sph[n] = sin(PHI2[n]);
                    m_w[n] = AREA2[n];
                }
                break;
            case 3:
                m_nint=NST3;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA3[n]);
                    m_sth[n] = sin(THETA3[n]);
                    m_cph[n] = cos(PHI3[n]);
                    m_sph[n] = sin(PHI3[n]);
                    m_w[n] = AREA3[n];
                }
                break;
            case 4:
                m_nint=NST4;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA4[n]);
                    m_sth[n] = sin(THETA4[n]);
                    m_cph[n] = cos(PHI4[n]);
                    m_sph[n] = sin(PHI4[n]);
                    m_w[n] = AREA4[n];
                }
                break;
            case 5:
                m_nint=NST5;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA5[n]);
                    m_sth[n] = sin(THETA5[n]);
                    m_cph[n] = cos(PHI5[n]);
                    m_sph[n] = sin(PHI5[n]);
                    m_w[n] = AREA5[n];
                }
                break;
            case 6:
                m_nint=NST6;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA6[n]);
                    m_sth[n] = sin(THETA6[n]);
                    m_cph[n] = cos(PHI6[n]);
                    m_sph[n] = sin(PHI6[n]);
                    m_w[n] = AREA6[n];
                }
                break;
            case 7:
                m_nint=NST7;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA7[n]);
                    m_sth[n] = sin(THETA7[n]);
                    m_cph[n] = cos(PHI7[n]);
                    m_sph[n] = sin(PHI7[n]);
                    m_w[n] = AREA7[n];
                }
                break;
            case 8:
                m_nint=NST8;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA8[n]);
                    m_sth[n] = sin(THETA8[n]);
                    m_cph[n] = cos(PHI8[n]);
                    m_sph[n] = sin(PHI8[n]);
                    m_w[n] = AREA8[n];
                }
                break;
            case 9:
                m_nint=NST9;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA9[n]);
                    m_sth[n] = sin(THETA9[n]);
                    m_cph[n] = cos(PHI9[n]);
                    m_sph[n] = sin(PHI9[n]);
                    m_w[n] = AREA9[n];
                }
                break;
                
            default:
                throw MaterialError("resolution must 1-9");
                break;
        }
        
		m_bfirst = false;
	}
    
    // also initialize the parent class
    FEFiberIntegrationScheme::Init();
}

//-----------------------------------------------------------------------------
mat3ds FEFiberIntegrationTriangle::Stress(FEMaterialPoint& mp)
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
		s += m_pFmat->Stress(mp)*(R*m_w[n]);
	}
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberIntegrationTriangle::Tangent(FEMaterialPoint& mp)
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
		c += m_pFmat->Tangent(mp)*(R*m_w[n]);
	}
	
	return c;
}

//-----------------------------------------------------------------------------
double FEFiberIntegrationTriangle::StrainEnergyDensity(FEMaterialPoint& mp)
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
		sed += m_pFmat->StrainEnergyDensity(mp)*(R*m_w[n]);
	}
	return sed;
}

//-----------------------------------------------------------------------------
void FEFiberIntegrationTriangle::IntegratedFiberDensity(double& IFD)
{
    // initialize integrated fiber density distribution
    IFD = 1;
    
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
}
