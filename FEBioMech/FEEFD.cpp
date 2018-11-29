//
//  FEEFD.cpp
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 5/18/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#include "stdafx.h"
#include "FEEFD.h"

// gaussian quadrature
/*const int nint = 2;
const double gp[nint] = {
    -sqrt(3.)/3.,
    sqrt(3.)/3.
};
const double gw[nint] = {
    1,
    1
};*/
/*const int nint = 3;
const double gp[nint] = {
    -sqrt(3./5.),
    0.,
    sqrt(3./5.)
};
const double gw[nint] = {
    5./9.,
    8./9.,
    5./9.
};*/
/*const int nint = 4;
const double gp[nint] = {
    -sqrt((3+2*sqrt(6./5.))/7.),
    -sqrt((3-2*sqrt(6./5.))/7.),
    sqrt((3-2*sqrt(6./5.))/7.),
    sqrt((3+2*sqrt(6./5.))/7.)
};
const double gw[nint] = {
    (18-sqrt(30.))/36.,
    (18+sqrt(30.))/36.,
    (18+sqrt(30.))/36.,
    (18-sqrt(30.))/36.
};*/

/*
const int nint = 5;
const double gp[nint] = {
    -1./3.*sqrt(5+2*sqrt(10./7.)),
    -1./3.*sqrt(5-2*sqrt(10./7.)),
    0.,
    1./3.*sqrt(5-2*sqrt(10./7.)),
    1./3.*sqrt(5+2*sqrt(10./7.))
};
const double gw[nint] = {
    (322-13*sqrt(70.))/900.,
    (322+13*sqrt(70.))/900.,
    128./225.,
    (322+13*sqrt(70.))/900.,
    (322-13*sqrt(70.))/900.
};

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FEEllipsoidalFiberDistribution
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEEFD, FEElasticMaterial)
	ADD_PARAMETER(m_beta, 3, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETER(m_ksi , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
mat3ds FEEFD::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
    
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
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);
	mat3d QT = Q.transpose();
    
    vec3d n0e, n0a, nt;
    double In, Wl, wn;
    
    // tension along all three eigenvectors (all directions)
    if (lC[0] > 1+eps) {
        ksia = 0;
        ksib = 1;
        dksi = (ksib - ksia)/2;
        sksi = (ksib + ksia)/2;
        
        // loop over all integration points
        for (int i=0; i<2*nint; ++i) {
            theta = i*pi/nint;
            for (int j=0; j<nint; ++j) {
                ksi = sksi + dksi*gp[j];
                wn = gw[j]*pi/nint*dksi;
                phi = acos(ksi);

                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);

                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                
                // calculate material coefficients
                double ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
                double beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));
                
                // get the global spatial fiber direction in current configuration
                nt = F*n0e;
                
                // Calculate In = n0e*C*n0e
                In = nt*nt;
                
                // calculate strain energy derivative
                Wl = beta*ksi*pow(In - 1.0, beta-1.0);
                
                // calculate the stress
                s += dyad(nt)*(Wl*wn);
            }
        }
    }
    // tension along one eigenvector and compression along other two
    else if (lC[1] <= 1+eps) {
        // loop over all integration points
        for (int i=0; i<2*nint; ++i) {
            theta = i*pi/nint;
            ksia = sqrt(1-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)))/sqrt(lC[2]-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)));
            ksib = 1;
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
            for (int j=0; j<nint; ++j) {
                ksi = sksi + dksi*gp[j];
                wn = gw[j]*pi/nint*dksi;
                phi = acos(ksi);
               
                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);
                
                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                
                // calculate material coefficients
                double ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
                double beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));
                
                // get the global spatial fiber direction in current configuration
                nt = F*n0e;
                
                // Calculate In = n0e*C*n0e
                In = nt*nt;
                
                // calculate strain energy derivative
                Wl = beta*ksi*pow(In - 1.0, beta-1.0);
                
                // calculate the stress
                s += dyad(nt)*(Wl*wn);
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
        for (int i=0; i<2*nint; ++i) {
            theta = i*pi/nint;
            ksia = 0;
            ksib = sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-1)/sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-lC[2]);
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
            for (int j=0; j<nint; ++j) {
                ksi = sksi + dksi*gp[j];
                wn = gw[j]*pi/nint*dksi;
                phi = acos(ksi);
                
                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);
                
                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                
                // calculate material coefficients
                double ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
                double beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));
                
                // get the global spatial fiber direction in current configuration
                nt = F*n0e;
                
                // Calculate In = n0e*C*n0e
                In = nt*nt;
                
                // calculate strain energy derivative
                Wl = beta*ksi*pow(In - 1.0, beta-1.0);
                
                // calculate the stress
                s += dyad(nt)*(Wl*wn);
            }
        }
    }
        
	// we multiply by two to add contribution from other half-sphere
	return s*(4.0/J);
}

//-----------------------------------------------------------------------------
tens4ds FEEFD::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
    
    // right Cauchy-Green tensor and its eigenvalues & eigenvectors
    mat3ds C = pt.RightCauchyGreen();
    double lC[3];
    vec3d vC[3];
    C.eigen(lC, vC);
    
    // initialize tangent tensor
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	tens4ds N4;
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
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);
	mat3d QT = Q.transpose();
    
    vec3d n0e, n0a, nt;
    double In, Wll, wn;
    
    // tension along all three eigenvectors (all directions)
    if (lC[0] > 1+eps) {
        ksia = 0;
        ksib = 1;
        dksi = (ksib - ksia)/2;
        sksi = (ksib + ksia)/2;
        
        // loop over all integration points
        for (int i=0; i<2*nint; ++i) {
            theta = i*pi/nint;
            for (int j=0; j<nint; ++j) {
                ksi = sksi + dksi*gp[j];
                wn = gw[j]*pi/nint*dksi;
                phi = acos(ksi);
                
                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);
                
                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                
                // calculate material coefficients
                double ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
                double beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));
                
                // get the global spatial fiber direction in current configuration
                nt = F*n0e;
                
                // Calculate In = n0e*C*n0e
                In = nt*nt;
                
                // calculate strain energy derivative
                Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
                
                N2 = dyad(nt);
                N4 = dyad1s(N2);
                
                c += N4*(Wll*wn);
            }
        }
    }
    // tension along one eigenvector and compression along other two
    else if (lC[1] <= 1+eps) {
        // loop over all integration points
        for (int i=0; i<2*nint; ++i) {
            theta = i*pi/nint;
            ksia = sqrt(1-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)))/sqrt(lC[2]-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)));
            ksib = 1;
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
            for (int j=0; j<nint; ++j) {
                ksi = sksi + dksi*gp[j];
                wn = gw[j]*pi/nint*dksi;
                phi = acos(ksi);
                
                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);
                
                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                
                // calculate material coefficients
                double ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
                double beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));
                
                // get the global spatial fiber direction in current configuration
                nt = F*n0e;
                
                // Calculate In = n0e*C*n0e
                In = nt*nt;
                
                // calculate strain energy derivative
                Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
                
                N2 = dyad(nt);
                N4 = dyad1s(N2);
                
                c += N4*(Wll*wn);
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
        for (int i=0; i<2*nint; ++i) {
            theta = i*pi/nint;
            ksia = 0;
            ksib = sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-1)/sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-lC[2]);
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
            for (int j=0; j<nint; ++j) {
                ksi = sksi + dksi*gp[j];
                wn = gw[j]*pi/nint*dksi;
                phi = acos(ksi);
                
                // set fiber direction in global coordinate system
                n0e = vC[0]*(cos(theta)*sin(phi)) + vC[1]*(sin(theta)*sin(phi)) + vC[2]*cos(phi);
                
                // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
                n0a = QT*n0e;
                
                // calculate material coefficients
                double ksi  = 1.0 / sqrt(SQR(n0a.x / m_ksi [0]) + SQR(n0a.y / m_ksi [1]) + SQR(n0a.z / m_ksi [2]));
                double beta = 1.0 / sqrt(SQR(n0a.x / m_beta[0]) + SQR(n0a.y / m_beta[1]) + SQR(n0a.z / m_beta[2]));
                
                // get the global spatial fiber direction in current configuration
                nt = F*n0e;
                
                // Calculate In = n0e*C*n0e
                In = nt*nt;
                
                // calculate strain energy derivative
                Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);
                
                N2 = dyad(nt);
                N4 = dyad1s(N2);
                
                c += N4*(Wll*wn);
            }
        }
    }
    
	// we multiply by two to add contribution from other half-sphere
	return c*(2.0*4.0/J);
}
*/
