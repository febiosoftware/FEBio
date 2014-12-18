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
                
            case 20:
                m_nint=NST20;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA20[n]);
                    m_sth[n] = sin(THETA20[n]);
                    m_cph[n] = cos(PHI20[n]);
                    m_sph[n] = sin(PHI20[n]);
                    m_w[n] = AREA20[n];
                }
                break;
            case 34:
                m_nint=NST34;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA34[n]);
                    m_sth[n] = sin(THETA34[n]);
                    m_cph[n] = cos(PHI34[n]);
                    m_sph[n] = sin(PHI34[n]);
                    m_w[n] = AREA34[n];
                }
                break;
            case 60:
                m_nint=NST60;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA60[n]);
                    m_sth[n] = sin(THETA60[n]);
                    m_cph[n] = cos(PHI60[n]);
                    m_sph[n] = sin(PHI60[n]);
                    m_w[n] = AREA60[n];
                }
                break;
            case 74:
                m_nint=NST74;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA74[n]);
                    m_sth[n] = sin(THETA74[n]);
                    m_cph[n] = cos(PHI74[n]);
                    m_sph[n] = sin(PHI74[n]);
                    m_w[n] = AREA74[n];
                }
                break;
            case 196:
                m_nint=NST196;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA196[n]);
                    m_sth[n] = sin(THETA196[n]);
                    m_cph[n] = cos(PHI196[n]);
                    m_sph[n] = sin(PHI196[n]);
                    m_w[n] = AREA196[n];
                }
                break;
            case 210:
                m_nint=NST210;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA210[n]);
                    m_sth[n] = sin(THETA210[n]);
                    m_cph[n] = cos(PHI210[n]);
                    m_sph[n] = sin(PHI210[n]);
                    m_w[n] = AREA210[n];
                }
                break;
                
            case 396:
                m_nint=NST396;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA396[n]);
                    m_sth[n] = sin(THETA396[n]);
                    m_cph[n] = cos(PHI396[n]);
                    m_sph[n] = sin(PHI396[n]);
                    m_w[n] = AREA396[n];
                }
                break;
                
            case 410:
                m_nint=NST410;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA410[n]);
                    m_sth[n] = sin(THETA410[n]);
                    m_cph[n] = cos(PHI410[n]);
                    m_sph[n] = sin(PHI410[n]);
                    m_w[n] = AREA410[n];
                }
                break;
            case 596:
                m_nint=NST596;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA596[n]);
                    m_sth[n] = sin(THETA596[n]);
                    m_cph[n] = cos(PHI596[n]);
                    m_sph[n] = sin(PHI596[n]);
                    m_w[n] = AREA596[n];
                }
                break;
            case 610:
                m_nint=NST610;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA610[n]);
                    m_sth[n] = sin(THETA610[n]);
                    m_cph[n] = cos(PHI610[n]);
                    m_sph[n] = sin(PHI610[n]);
                    m_w[n] = AREA610[n];
                }
                break;
            case 796:
                m_nint=NST796;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA796[n]);
                    m_sth[n] = sin(THETA796[n]);
                    m_cph[n] = cos(PHI796[n]);
                    m_sph[n] = sin(PHI796[n]);
                    m_w[n] = AREA796[n];
                }
                break;
            case 810:
                m_nint=NST810;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA810[n]);
                    m_sth[n] = sin(THETA810[n]);
                    m_cph[n] = cos(PHI810[n]);
                    m_sph[n] = sin(PHI810[n]);
                    m_w[n] = AREA810[n];
                }
                break;
                
            case 996:
                m_nint=NST996;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA996[n]);
                    m_sth[n] = sin(THETA996[n]);
                    m_cph[n] = cos(PHI996[n]);
                    m_sph[n] = sin(PHI996[n]);
                    m_w[n] = AREA996[n];
                }
                break;
                
            case 1010:
                m_nint=NST1010;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA1010[n]);
                    m_sth[n] = sin(THETA1010[n]);
                    m_cph[n] = cos(PHI1010[n]);
                    m_sph[n] = sin(PHI1010[n]);
                    m_w[n] = AREA1010[n];
                }
                break;
                
            case 1196:
                m_nint=NST1196;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA1196[n]);
                    m_sth[n] = sin(THETA1196[n]);
                    m_cph[n] = cos(PHI1196[n]);
                    m_sph[n] = sin(PHI1196[n]);
                    m_w[n] = AREA1196[n];
                }
                break;
                
            case 1210:
                m_nint=NST1210;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA1210[n]);
                    m_sth[n] = sin(THETA1210[n]);
                    m_cph[n] = cos(PHI1210[n]);
                    m_sph[n] = sin(PHI1210[n]);
                    m_w[n] = AREA1210[n];
                }
                break;
            case 1396:
                m_nint=NST1396;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA1396[n]);
                    m_sth[n] = sin(THETA1396[n]);
                    m_cph[n] = cos(PHI1396[n]);
                    m_sph[n] = sin(PHI1396[n]);
                    m_w[n] = AREA1396[n];
                }
                break;
            case 1410:
                m_nint=NST1410;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA1410[n]);
                    m_sth[n] = sin(THETA1410[n]);
                    m_cph[n] = cos(PHI1410[n]);
                    m_sph[n] = sin(PHI1410[n]);
                    m_w[n] = AREA1410[n];
                }
                break;
            case 1596:
                m_nint=NST1596;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA1596[n]);
                    m_sth[n] = sin(THETA1596[n]);
                    m_cph[n] = cos(PHI1596[n]);
                    m_sph[n] = sin(PHI1596[n]);
                    m_w[n] = AREA1596[n];
                }
                break;
            case 1610:
                m_nint=NST1610;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA1610[n]);
                    m_sth[n] = sin(THETA1610[n]);
                    m_cph[n] = cos(PHI1610[n]);
                    m_sph[n] = sin(PHI1610[n]);
                    m_w[n] = AREA1610[n];
                }
                break;
            case 1796:
                m_nint=NST1796;
                for (int n=0; n<m_nint; ++n)
                {
                    m_cth[n] = cos(THETA1796[n]);
                    m_sth[n] = sin(THETA1796[n]);
                    m_cph[n] = cos(PHI1796[n]);
                    m_sph[n] = sin(PHI1796[n]);
                    m_w[n] = AREA1796[n];
                }
                break;
                
                //            default:
                //                throw MaterialError("resolutio must 1-9");
                //                break;
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
