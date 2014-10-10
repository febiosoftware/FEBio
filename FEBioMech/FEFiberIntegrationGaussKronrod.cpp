//
//  FEFiberIntegrationGaussKronrod.cpp
//  FEBioXCode4
//
//

#include "stdafx.h"
#include "FEFiberIntegrationGaussKronrod.h"
#include "gausskronrod.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// FEFiberIntegrationGaussKronrod
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberIntegrationGaussKronrod, FEFiberIntegrationScheme)
    ADD_PARAMETER(m_nph, FE_PARAM_INT, "nph");
    ADD_PARAMETER(m_nth, FE_PARAM_INT, "nth");
END_PARAMETER_LIST();

void FEFiberIntegrationGaussKronrod::Init()
{
	if (m_nph < 1) throw MaterialError("nph must be strictly greater than zero.");
	if (m_nth < 1) throw MaterialError("nth must be strictly greater than zero.");
    
	m_bfirst = true;
	
	if (m_bfirst)
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
		m_bfirst = false;
    }
    
    // also initialize the parent class
    FEFiberIntegrationScheme::Init();
}

//-----------------------------------------------------------------------------
mat3ds FEFiberIntegrationGaussKronrod::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the element's local coordinate system
	mat3d QT = (pt.m_Q).transpose();
	
    // Lagrangian strain tensor and its eigenvalues & eigenvectors
    mat3ds E = pt.Strain();
    mat3ds C = pt.RightCauchyGreen();
    double lE[3], lC[3];
    vec3d vE[3], vC[3];
    E.eigen2(lE,vE);
    C.eigen2(lC,vC);
    
    // initialize stress tensor
	mat3ds s;
	s.zero();
    
    // check if there is no tension (C case)
    if (lE[2] <= 0) {
        return s;
    }
    
    // check remaining stretch states
    double phi, theta;
    double ksia, ksib, dksi, sksi, ksi;
    double pi = 4*atan(1.0);
    double dth = 2*pi/m_nth;
    
    vec3d n0e, n0a;
    double wn;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        double den = lE[0]*SQR(cos(theta)) + lE[1]*SQR(sin(theta));
        
        // check for T, TC or CT cases
        
        // T case
        if (lE[0] >= 0) {
            ksia = 0;
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
        // TC case
        else if (lE[1] <= 0) {
            ksia = sqrt(den/(den-lE[2]));
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
        // CT case
        else {
            ksia = 0;
            ksib = sqrt(den/(den-lE[2]));
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
tens4ds FEFiberIntegrationGaussKronrod::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the element's local coordinate system
    mat3d QT = (pt.m_Q).transpose();
    
    // Lagrangian strain tensor and its eigenvalues & eigenvectors
    mat3ds E = pt.Strain();
    mat3ds C = pt.RightCauchyGreen();
    double lE[3], lC[3];
    vec3d vE[3], vC[3];
    E.eigen2(lE,vE);
    C.eigen2(lC,vC);
    
    // initialize tangent tensor
    tens4ds c;
    c.zero();
    
    // check if there is no tension (C case)
    if (lE[2] <= 0) {
        return c;
    }
    
    // check remaining stretch states
    double phi, theta;
    double ksia, ksib, dksi, sksi, ksi;
    double pi = 4*atan(1.0);
    double dth = 2*pi/m_nth;
    
    vec3d n0e, n0a;
    double wn;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        double den = lE[0]*SQR(cos(theta)) + lE[1]*SQR(sin(theta));
        
        // check for T, TC or CT cases
        
        // T case
        if (lE[0] >= 0) {
            ksia = 0;
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
        // TC case
        else if (lE[1] <= 0) {
            ksia = sqrt(den/(den-lE[2]));
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
        // CT case
        else {
            ksia = 0;
            ksib = sqrt(den/(den-lE[2]));
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
double FEFiberIntegrationGaussKronrod::StrainEnergyDensity(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the element's local coordinate system
    mat3d QT = (pt.m_Q).transpose();
    
    // Lagrangian strain tensor and its eigenvalues & eigenvectors
    mat3ds E = pt.Strain();
    mat3ds C = pt.RightCauchyGreen();
    double lE[3], lC[3];
    vec3d vE[3], vC[3];
    E.eigen2(lE,vE);
    C.eigen2(lC,vC);
    
    // initialize strain energy density
    double sed = 0.0;
    
    // check if there is no tension (C case)
    if (lE[2] <= 0) {
        return sed;
    }
    
    // check remaining stretch states
    double phi, theta;
    double ksia, ksib, dksi, sksi, ksi;
    double pi = 4*atan(1.0);
    double dth = 2*pi/m_nth;
    
    vec3d n0e, n0a;
    double wn;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        double den = lE[0]*SQR(cos(theta)) + lE[1]*SQR(sin(theta));
        
        // check for T, TC or CT cases
        
        // T case
        if (lE[0] >= 0) {
            ksia = 0;
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
        // TC case
        else if (lE[1] <= 0) {
            ksia = sqrt(den/(den-lE[2]));
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
        // CT case
        else {
            ksia = 0;
            ksib = sqrt(den/(den-lE[2]));
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
void FEFiberIntegrationGaussKronrod::IntegratedFiberDensity(double& IFD)
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
