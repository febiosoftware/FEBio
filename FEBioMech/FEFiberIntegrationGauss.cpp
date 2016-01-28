//
//  FEFiberIntegrationGauss.cpp
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 11/19/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#include "FEFiberIntegrationGauss.h"
#include "gauss.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FEFiberIntegrationGauss
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberIntegrationGauss, FEFiberIntegrationScheme)
	ADD_PARAMETER2(m_nph, FE_PARAM_INT, FE_RANGE_GREATER(0), "nph");
	ADD_PARAMETER2(m_nth, FE_PARAM_INT, FE_RANGE_GREATER(0), "nth");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEFiberIntegrationGauss::Serialize(DumpStream& ar)
{
	FEFiberIntegrationScheme::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_gp << m_gw;
	}
	else
	{
		ar >> m_gp >> m_gw;
	}
}

//-----------------------------------------------------------------------------
bool FEFiberIntegrationGauss::Init()
{
    switch (m_nph) {
        case 1:
            m_gp.assign(gp1, gp1+nint1);
            m_gw.assign(gw1, gw1+nint1);
            break;
        case 2:
            m_gp.assign(gp2, gp2+nint2);
            m_gw.assign(gw2, gw2+nint2);
            break;
        case 3:
            m_gp.assign(gp3, gp3+nint3);
            m_gw.assign(gw3, gw3+nint3);
            break;
        case 4:
            m_gp.assign(gp4, gp4+nint4);
            m_gw.assign(gw4, gw4+nint4);
            break;
        case 5:
            m_gp.assign(gp5, gp5+nint5);
            m_gw.assign(gw5, gw5+nint5);
            break;
        case 6:
            m_gp.assign(gp6, gp6+nint6);
            m_gw.assign(gw6, gw6+nint6);
            break;
        case 7:
            m_gp.assign(gp7, gp7+nint7);
            m_gw.assign(gw7, gw7+nint7);
            break;
        case 8:
            m_gp.assign(gp8, gp8+nint8);
            m_gw.assign(gw8, gw8+nint8);
            break;
        case 9:
            m_gp.assign(gp9, gp9+nint9);
            m_gw.assign(gw9, gw9+nint9);
            break;
        case 10:
            m_gp.assign(gp10, gp10+nint10);
            m_gw.assign(gw10, gw10+nint10);
            break;
        default:
            return MaterialError("nint must not exceed 10.");
            break;
    }
    
    // also initialize the parent class
    return FEFiberIntegrationScheme::Init();
}

//-----------------------------------------------------------------------------
mat3ds FEFiberIntegrationGauss::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();
	
    // right Cauchy-Green tensor and its eigenvalues & eigenvectors
    mat3ds C = pt.RightCauchyGreen();
    double lC[3];
    vec3d vC[3];
    C.eigen(lC, vC);
    
    // initialize stress tensor
	mat3ds s;
	s.zero();
    
    // check if there is no tension
    const double eps = 1e-9;
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
                s += m_pFmat->Stress(mp)*(R*wn);
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
                s += m_pFmat->Stress(mp)*(R*wn);
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
                s += m_pFmat->Stress(mp)*(R*wn);
            }
        }
    }
    
	// we multiply by two to add contribution from other half-sphere
	return s*(2.0);
}

//-----------------------------------------------------------------------------
tens4ds FEFiberIntegrationGauss::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();
	
    // right Cauchy-Green tensor and its eigenvalues & eigenvectors
    mat3ds C = pt.RightCauchyGreen();
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
                
                // calculate the stress
                c += m_pFmat->Tangent(mp)*(R*wn);
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
                c += m_pFmat->Tangent(mp)*(R*wn);
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
                c += m_pFmat->Tangent(mp)*(R*wn);
            }
        }
    }
    
	// we multiply by two to add contribution from other half-sphere
	return c*(2.0);
}

//-----------------------------------------------------------------------------
double FEFiberIntegrationGauss::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();
	
    // right Cauchy-Green tensor and its eigenvalues & eigenvectors
    mat3ds C = pt.RightCauchyGreen();
    double lC[3];
    vec3d vC[3];
    C.eigen(lC, vC);
    
    // initialize strain energy density
    double sed = 0.0;
    
    // check if there is no tension
    const double eps = 1e-9;
    if ((lC[0] <= 1+eps) && (lC[1] <= 1+eps) && (lC[2] <= 1+eps)) {
        return sed;
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
                sed += m_pFmat->StrainEnergyDensity(mp)*(R*wn);
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
                sed += m_pFmat->StrainEnergyDensity(mp)*(R*wn);
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
                sed += m_pFmat->StrainEnergyDensity(mp)*(R*wn);
            }
        }
    }
    
	// we multiply by two to add contribution from other half-sphere
	return sed*2.0;
}

//-----------------------------------------------------------------------------
void FEFiberIntegrationGauss::IntegratedFiberDensity(double& IFD)
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
