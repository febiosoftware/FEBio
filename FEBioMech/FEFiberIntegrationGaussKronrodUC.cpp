//
//  FEFiberIntegrationGaussKronrodUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEFiberIntegrationGaussKronrodUC.h"
#include "FEContinuousFiberDistributionUC.h"
#include "gausskronrod.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FEFiberIntegrationGaussKronrodUC
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberIntegrationGaussKronrodUC, FEFiberIntegrationSchemeUC)
    ADD_PARAMETER(m_nph, FE_PARAM_INT, "nph");
    ADD_PARAMETER(m_nth, FE_PARAM_INT, "nth");
END_PARAMETER_LIST();

void FEFiberIntegrationGaussKronrodUC::Init()
{
	if (m_nph < 1) throw MaterialError("nph must be strictly greater than zero.");
	if (m_nth < 1) throw MaterialError("nth must be strictly greater than zero.");
    
	static bool bfirst = true;
	
	if (bfirst)
	{
        switch (m_nph) {
            case 7:
                m_gp.assign(gp7, gp7+nint7);
                m_gw.assign(gw7, gw7+nint7);
                break;
            case 11:
                m_gp.assign(gp11, gp11+nint11);
                m_gw.assign(gw11, gw11+nint11);
                break;
            case 15:
                m_gp.assign(gp15, gp15+nint15);
                m_gw.assign(gw15, gw15+nint15);
                break;
            case 19:
                m_gp.assign(gp19, gp19+nint19);
                m_gw.assign(gw19, gw19+nint19);
                break;
            case 23:
                m_gp.assign(gp23, gp23+nint23);
                m_gw.assign(gw23, gw23+nint23);
                break;
            case 27:
                m_gp.assign(gp27, gp27+nint27);
                m_gw.assign(gw27, gw27+nint27);
                break;
            default:
                throw MaterialError("nph must 7,11,15,19,23,27.");
                break;
        }
		bfirst = false;
    }
    
    // also initialize the parent class
    FEFiberIntegrationSchemeUC::Init();
}

//-----------------------------------------------------------------------------
mat3ds FEFiberIntegrationGaussKronrodUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEContinuousFiberDistributionUC* pcfd = dynamic_cast<FEContinuousFiberDistributionUC*>(GetParent());
	
	// get the element's local coordinate system
	mat3d QT = (pcfd->LocalMatAxes()*pt.m_Q).transpose();
	
    // right Cauchy-Green tensor and its eigenvalues & eigenvectors
    mat3ds C = pt.DevRightCauchyGreen();
    double lC[3];
    vec3d vC[3];
    C.eigen(lC, vC);
    
    // initialize stress tensor
	mat3ds s;
	s.zero();
    
    // check if there is no tension
    const double eps = 1.e-9;
    if ((lC[0] <= 1+eps) && (lC[1] <= 1+eps) && (lC[2] <= 1+eps)) {
        return s;
    }
    
    // bubble sort eigenvalues & eigenvectors from smallest to largest
    double ltmp;
    vec3d vtmp;
    bool swp = true;
    while (swp) {
        swp = false;
        for (int i=1; i<3; ++i) {
            int j = i-1;
            if (lC[i] < lC[j]) {
                ltmp = lC[i]; vtmp = vC[i];
                lC[i] = lC[j]; vC[i] = vC[j];
                lC[j] = ltmp; vC[j] = vtmp;
                swp = true;
            }
        }
    }
    
    // check remaining stretch states
    double phi, theta;
    double ksia, ksib, dksi, sksi, ksi;
    double pi = 4*atan(1.0);
    double dth = 2*pi/m_nth;
    
    vec3d n0e, n0a;
    double wn;
    
    // tension along all three eigenvectors (all directions)
    if (lC[0] > 1+eps) {
        ksia = 0;
        ksib = 1;
        dksi = (ksib - ksia)/2;
        sksi = (ksib + ksia)/2;
        
        // loop over all integration points
        for (int i=0; i<m_nth; ++i) {
            theta = i*dth;
            for (int j=0; j<m_nph; ++j) {
                ksi = sksi + dksi*m_gp[j];
                wn = m_gw[j]*dth*dksi;
                phi = acos(ksi);
                
                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);
                m_pFmat->SetFiberDirection(mp, n0e);
                
                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                double R = m_pFDD->FiberDensity(n0a);
                
                // calculate the stress
                s += m_pFmat->DevStress(mp)*(R*wn);
            }
        }
    }
    // tension along one eigenvector and compression along other two
    else if (lC[1] <= 1+eps) {
        // loop over all integration points
        for (int i=0; i<m_nth; ++i) {
            theta = i*dth;
            ksia = sqrt(1-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)))/sqrt(lC[2]-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)));
            ksib = 1;
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
            for (int j=0; j<m_nph; ++j) {
                ksi = sksi + dksi*m_gp[j];
                wn = m_gw[j]*dth*dksi;
                phi = acos(ksi);
                
                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);
                m_pFmat->SetFiberDirection(mp, n0e);
                
                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                double R = m_pFDD->FiberDensity(n0a);
                
                // calculate the stress
                s += m_pFmat->DevStress(mp)*(R*wn);
            }
        }
    }
    // tension along two eigenvectors and compression along third
    else {
        // swap first and last eigenvalues/eigenvectors to maintain consistency in formulas
        ltmp = lC[2]; vtmp = vC[2];
        lC[2] = lC[0]; vC[2] = vC[0];
        lC[0] = ltmp; vC[0] = vtmp;
        
        // loop over all integration points
        for (int i=0; i<m_nth; ++i) {
            theta = i*dth;
            ksia = 0;
            ksib = sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-1)/sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-lC[2]);
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
            for (int j=0; j<m_nph; ++j) {
                ksi = sksi + dksi*m_gp[j];
                wn = m_gw[j]*dth*dksi;
                phi = acos(ksi);
                
                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);
                m_pFmat->SetFiberDirection(mp, n0e);
                
                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                double R = m_pFDD->FiberDensity(n0a);
                
                // calculate the stress
                s += m_pFmat->DevStress(mp)*(R*wn);
            }
        }
    }
    
	// we multiply by two to add contribution from other half-sphere
    // we don't need to evaluate the deviatoric part since we already summed up deviatoric fiber stresses
	return s*(2.0);
}

//-----------------------------------------------------------------------------
tens4ds FEFiberIntegrationGaussKronrodUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEContinuousFiberDistributionUC* pcfd = dynamic_cast<FEContinuousFiberDistributionUC*>(GetParent());
	
	// get the element's local coordinate system
	mat3d QT = (pcfd->LocalMatAxes()*pt.m_Q).transpose();
	
    // right Cauchy-Green tensor and its eigenvalues & eigenvectors
    mat3ds C = pt.DevRightCauchyGreen();
    double lC[3];
    vec3d vC[3];
    C.eigen(lC, vC);
    
    // initialize stress tensor
	tens4ds c;
	c.zero();
    
    // check if there is no tension
    const double eps = 1e-9;
    if ((lC[0] <= 1+eps) && (lC[1] <= 1+eps) && (lC[2] <= 1+eps)) {
        return c;
    }
    
    // bubble sort eigenvalues & eigenvectors from smallest to largest
    double ltmp;
    vec3d vtmp;
    bool swp = true;
    while (swp) {
        swp = false;
        for (int i=1; i<3; ++i) {
            int j = i-1;
            if (lC[i] < lC[j]) {
                ltmp = lC[i]; vtmp = vC[i];
                lC[i] = lC[j]; vC[i] = vC[j];
                lC[j] = ltmp; vC[j] = vtmp;
                swp = true;
            }
        }
    }
    
    // check remaining stretch states
    double phi, theta;
    double ksia, ksib, dksi, sksi, ksi;
    double pi = 4*atan(1.0);
    double dth = 2*pi/m_nth;
    
    vec3d n0e, n0a;
    double wn;
    
    // tension along all three eigenvectors (all directions)
    if (lC[0] > 1+eps) {
        ksia = 0;
        ksib = 1;
        dksi = (ksib - ksia)/2;
        sksi = (ksib + ksia)/2;
        
        // loop over all integration points
        for (int i=0; i<m_nth; ++i) {
            theta = i*dth;
            for (int j=0; j<m_nph; ++j) {
                ksi = sksi + dksi*m_gp[j];
                wn = m_gw[j]*dth*dksi;
                phi = acos(ksi);
                
                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);
                m_pFmat->SetFiberDirection(mp, n0e);
                
                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                double R = m_pFDD->FiberDensity(n0a);
                
                // calculate the tangent
                c += m_pFmat->DevTangent(mp)*(R*wn);
            }
        }
    }
    // tension along one eigenvector and compression along other two
    else if (lC[1] <= 1+eps) {
        // loop over all integration points
        for (int i=0; i<m_nth; ++i) {
            theta = i*dth;
            ksia = sqrt(1-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)))/sqrt(lC[2]-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)));
            ksib = 1;
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
            for (int j=0; j<m_nph; ++j) {
                ksi = sksi + dksi*m_gp[j];
                wn = m_gw[j]*dth*dksi;
                phi = acos(ksi);
                
                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);
                m_pFmat->SetFiberDirection(mp, n0e);
                
                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                double R = m_pFDD->FiberDensity(n0a);
                
                // calculate the tangent
                c += m_pFmat->DevTangent(mp)*(R*wn);
            }
        }
    }
    // tension along two eigenvectors and compression along third
    else {
        // swap first and last eigenvalues/eigenvectors to maintain consistency in formulas
        ltmp = lC[2]; vtmp = vC[2];
        lC[2] = lC[0]; vC[2] = vC[0];
        lC[0] = ltmp; vC[0] = vtmp;
        
        // loop over all integration points
        for (int i=0; i<m_nth; ++i) {
            theta = i*dth;
            ksia = 0;
            ksib = sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-1)/sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-lC[2]);
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
            for (int j=0; j<m_nph; ++j) {
                ksi = sksi + dksi*m_gp[j];
                wn = m_gw[j]*dth*dksi;
                phi = acos(ksi);
                
                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);
                m_pFmat->SetFiberDirection(mp, n0e);
                
                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                double R = m_pFDD->FiberDensity(n0a);
                
                // calculate the tangent
                c += m_pFmat->DevTangent(mp)*(R*wn);
            }
        }
    }
    
	// we multiply by two to add contribution from other half-sphere
    // we don't need to evaluate the deviatoric part since we already summed up deviatoric fiber tangents
	return c*(2.0);
}

//-----------------------------------------------------------------------------
void FEFiberIntegrationGaussKronrodUC::IntegratedFiberDensity(double& IFD)
{
    // establish local basis
    vec3d a[3], n0a;
    a[0] = vec3d(1,0,0);
    a[1] = vec3d(0,1,0);
    a[2] = vec3d(0,0,1);
    
    // initialize integrated fiber density
    IFD = 1;
    double C = 0;
    
    double phi, theta;
    double ksia, ksib, dksi, sksi, ksi;
    double wn;
    double pi = 4*atan(1.0);
    double dth = 2*pi/m_nth;
    
    ksia = 0;
    ksib = 1;
    dksi = (ksib - ksia)/2;
    sksi = (ksib + ksia)/2;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        for (int j=0; j<m_nph; ++j) {
            ksi = sksi + dksi*m_gp[j];
            wn = m_gw[j]*dth*dksi;
            phi = acos(ksi);
            
            // set fiber direction in local coordinate system
            n0a = a[0]*(cos(theta)*sin(phi)) + a[1]*(sin(theta)*sin(phi)) + a[2]*cos(phi);
            // evaluate fiber density along this direction
            double R = m_pFDD->FiberDensity(n0a);
            
            // integrate fiber density
            C += R*wn;
        }
    }
    
	// we multiply by two to add contribution from other half-sphere
	IFD = C*2;
}
