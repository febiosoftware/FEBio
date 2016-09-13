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

class FEFiberIntegrationGauss::Iterator : public FEFiberIntegrationSchemeIterator
{
public:
	Iterator(FEMaterialPoint* mp, FEFiberIntegrationGauss::GRULE& rule)
	{
		m_ncase = -1;
		m_nth = rule.m_nth;
		m_nph = rule.m_nph;
		m_gp = rule.m_gp;
		m_gw = rule.m_gw;

		const double eps = 1e-9;
		if (mp)
		{
			FEElasticMaterialPoint& pt = *mp->ExtractData<FEElasticMaterialPoint>();
			// right Cauchy-Green tensor and its eigenvalues & eigenvectors
			// TODO: for uncoupled formulations we need to use the deviatoric version
			// mat3ds C = pt.DevRightCauchyGreen();
			mat3ds C = pt.RightCauchyGreen();
			C.eigen(lC, vC);

			// check if there is no tension
			if ((lC[0] <= 1 + eps) && (lC[1] <= 1 + eps) && (lC[2] <= 1 + eps)) {
				return;
			}

			// bubble sort eigenvalues & eigenvectors from smallest to largest
			double ltmp;
			vec3d vtmp;
			bool swp = true;
			while (swp) {
				swp = false;
				for (int i = 1; i<3; ++i) {
					int j = i - 1;
					if (lC[i] < lC[j]) {
						ltmp = lC[i]; vtmp = vC[i];
						lC[i] = lC[j]; vC[i] = vC[j];
						lC[j] = ltmp; vC[j] = vtmp;
						swp = true;
					}
				}
			}

			// check the other case
			if      (lC[0]  > 1 + eps) m_ncase = 0;
			else if (lC[1] <= 1 + eps) m_ncase = 1;
			else m_ncase = 2;
		}
		else
		{
			m_ncase = 0;
			vC[0] = vec3d(1,0,0);
			vC[1] = vec3d(0,1,0);
			vC[2] = vec3d(0,0,1);
		}

		// check remaining stretch states
		double pi = 4 * atan(1.0);
		dth = 2 * pi / m_nth;

		vec3d n0e, n0a;

		// tension along all three eigenvectors (all directions)
		if (m_ncase == 0) 
		{
			ksia = 0;
			ksib = 1;
			dksi = (ksib - ksia) / 2;
			sksi = (ksib + ksia) / 2;
		}
		// tension along one eigenvector and compression along other two
		else if (m_ncase == 1)
		{
			// nothing to do
			// ksia and ksib are update in FiberVector
		}
		// tension along two eigenvectors and compression along third
		else
		{
			// swap first and last eigenvalues/eigenvectors to maintain consistency in formulas
			double ltmp = lC[2]; vec3d vtmp = vC[2];
			lC[2] = lC[0]; vC[2] = vC[0];
			lC[0] = ltmp; vC[0] = vtmp;
		}

		i = 0; j = -1;
		i_old = -1;
		cth = 1.0; sth = 0.0;
		Next();
	}

	bool IsValid()
	{
		return (m_ncase != -1);
	}

	// move to the next integration point
	bool Next()
	{
		// make sure the iterator is valid
		if (m_ncase == -1) return false;

		// update loop counters
		j++;
		if (j >= m_nph)
		{
			j = 0;
			i++;
			if (i >= m_nth)
			{
				// all done
				m_ncase = -1;
				return false;
			}
		}

		// update vector
		if (i_old != i)
		{
			double theta = i*dth;
			cth = cos(theta);
			sth = sin(theta);

			if (m_ncase == 1)
			{
				ksia = sqrt(1 - lC[0] * SQR(cth) - lC[1] * SQR(sth)) / sqrt(lC[2] - lC[0] * SQR(cth) - lC[1] * SQR(sth));
				ksib = 1;
				dksi = (ksib - ksia) / 2;
				sksi = (ksib + ksia) / 2;
			}
			else if (m_ncase == 2)
			{
				ksia = 0;
				ksib = sqrt(lC[0] * SQR(cth) + lC[1] * SQR(sth) - 1) / sqrt(lC[0] * SQR(cth) + lC[1] * SQR(sth) - lC[2]);
				dksi = (ksib - ksia) / 2;
				sksi = (ksib + ksia) / 2;
			}
		}

		double ksi = sksi + dksi*m_gp[j];
		double sph = sqrt(1.0 - ksi*ksi); // = sin(acos(ksi));

		m_fiber = vC[0] * (cth*sph) + vC[1] * (sth*sph) + vC[2] * ksi;

		// we multiply by two to add contribution from other half-sphere
		m_weight = (m_gw[j] * dth*dksi)*2.0;

		i_old = i;

		return true;
	}

public:
	int m_ncase;
	double lC[3];
	vec3d vC[3];
	int	m_nth, m_nph;
	const double* m_gp;
	const double* m_gw;

	int i, j, i_old;
	double ksia, ksib, dksi, sksi;
	double dth;
	double cth, sth;
};

//-----------------------------------------------------------------------------
// FEFiberIntegrationGauss
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberIntegrationGauss, FEFiberIntegrationScheme)
	ADD_PARAMETER2(m_rule.m_nph, FE_PARAM_INT, FE_RANGE_GREATER(0), "nph");
	ADD_PARAMETER2(m_rule.m_nth, FE_PARAM_INT, FE_RANGE_GREATER(0), "nth");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEFiberIntegrationGauss::Serialize(DumpStream& ar)
{
	FEFiberIntegrationScheme::Serialize(ar);
	if ((ar.IsSaving() == false) && (ar.IsShallow() == false))
	{
		InitRule();
	}
}

//-----------------------------------------------------------------------------
bool FEFiberIntegrationGauss::InitRule()
{
	switch (m_rule.m_nph) {
        case 1:
			m_rule.m_gp = gp1;
			m_rule.m_gw = gw1;
            break;
        case 2:
			m_rule.m_gp = gp2;
			m_rule.m_gw = gw2;
            break;
        case 3:
			m_rule.m_gp = gp3;
			m_rule.m_gw = gw3;
            break;
        case 4:
			m_rule.m_gp = gp4;
			m_rule.m_gw = gw4;
            break;
        case 5:
			m_rule.m_gp = gp5;
			m_rule.m_gw = gw5;
            break;
        case 6:
			m_rule.m_gp = gp6;
			m_rule.m_gw = gw6;
            break;
        case 7:
			m_rule.m_gp = gp7;
			m_rule.m_gw = gw7;
            break;
        case 8:
			m_rule.m_gp = gp8;
			m_rule.m_gw = gw8;
            break;
        case 9:
			m_rule.m_gp = gp9;
			m_rule.m_gw = gw9;
            break;
        case 10:
			m_rule.m_gp = gp10;
			m_rule.m_gw = gw10;
            break;
        default:
            return false;
            break;
    }

	return true;
}

FEFiberIntegrationGauss::FEFiberIntegrationGauss(FEModel* pfem) : FEFiberIntegrationScheme(pfem)
{ 
	m_rule.m_nph = 5; 
	m_rule.m_nth = 2 * m_rule.m_nph;
}

FEFiberIntegrationGauss::~FEFiberIntegrationGauss()
{
}

bool FEFiberIntegrationGauss::Init()
{
	if (InitRule() == false) return MaterialError("nint must not exceed 10.");

    // also initialize the parent class
    return FEFiberIntegrationScheme::Init();
}

//-----------------------------------------------------------------------------
FEFiberIntegrationSchemeIterator* FEFiberIntegrationGauss::GetIterator(FEMaterialPoint* mp)
{
	return new Iterator(mp, m_rule);
}

/*
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
	double dth = 2 * pi / m_rule.m_nth;
    
    vec3d n0e, n0a;
    double wn;
    
    // tension along all three eigenvectors (all directions)
    if (lC[0] > 1+eps) {
        ksia = 0;
        ksib = 1;
        dksi = (ksib - ksia)/2;
        sksi = (ksib + ksia)/2;
        
        // loop over all integration points
		for (int i = 0; i<m_rule.m_nth; ++i) {
            theta = i*dth;
			for (int j = 0; j<m_rule.m_nph; ++j) {
				ksi = sksi + dksi*m_rule.m_gp[j];
				wn = m_rule.m_gw[j] * dth*dksi;
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
		for (int i = 0; i<m_rule.m_nth; ++i) {
            theta = i*dth;
            ksia = sqrt(1-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)))/sqrt(lC[2]-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)));
            ksib = 1;
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
			for (int j = 0; j<m_rule.m_nph; ++j) {
				ksi = sksi + dksi*m_rule.m_gp[j];
				wn = m_rule.m_gw[j] * dth*dksi;
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
		for (int i = 0; i<m_rule.m_nth; ++i) {
            theta = i*dth;
            ksia = 0;
            ksib = sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-1)/sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-lC[2]);
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
			for (int j = 0; j<m_rule.m_nph; ++j) {
				ksi = sksi + dksi*m_rule.m_gp[j];
				wn = m_rule.m_gw[j] * dth*dksi;
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
	double dth = 2 * pi / m_rule.m_nth;
    
    vec3d n0e, n0a;
    double wn;
    
    // tension along all three eigenvectors (all directions)
    if (lC[0] > 1+eps) {
        ksia = 0;
        ksib = 1;
        dksi = (ksib - ksia)/2;
        sksi = (ksib + ksia)/2;
        
        // loop over all integration points
		for (int i = 0; i<m_rule.m_nth; ++i) {
            theta = i*dth;
			for (int j = 0; j<m_rule.m_nph; ++j) {
				ksi = sksi + dksi*m_rule.m_gp[j];
				wn = m_rule.m_gw[j] * dth*dksi;
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
		for (int i = 0; i<m_rule.m_nth; ++i) {
            theta = i*dth;
            ksia = sqrt(1-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)))/sqrt(lC[2]-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)));
            ksib = 1;
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
			for (int j = 0; j<m_rule.m_nph; ++j) {
				ksi = sksi + dksi*m_rule.m_gp[j];
				wn = m_rule.m_gw[j] * dth*dksi;
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
		for (int i = 0; i<m_rule.m_nth; ++i) {
            theta = i*dth;
            ksia = 0;
            ksib = sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-1)/sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-lC[2]);
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
			for (int j = 0; j<m_rule.m_nph; ++j) {
				ksi = sksi + dksi*m_rule.m_gp[j];
				wn = m_rule.m_gw[j] * dth*dksi;
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
	double dth = 2 * pi / m_rule.m_nth;
    
    vec3d n0e, n0a;
    double wn;
    
    // tension along all three eigenvectors (all directions)
    if (lC[0] > 1+eps) {
        ksia = 0;
        ksib = 1;
        dksi = (ksib - ksia)/2;
        sksi = (ksib + ksia)/2;
        
        // loop over all integration points
		for (int i = 0; i<m_rule.m_nth; ++i) {
            theta = i*dth;
			for (int j = 0; j<m_rule.m_nph; ++j) {
				ksi = sksi + dksi*m_rule.m_gp[j];
				wn = m_rule.m_gw[j] * dth*dksi;
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
		for (int i = 0; i<m_rule.m_nth; ++i) {
            theta = i*dth;
            ksia = sqrt(1-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)))/sqrt(lC[2]-lC[0]*SQR(cos(theta))-lC[1]*SQR(sin(theta)));
            ksib = 1;
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
			for (int j = 0; j<m_rule.m_nph; ++j) {
				ksi = sksi + dksi*m_rule.m_gp[j];
				wn = m_rule.m_gw[j] * dth*dksi;
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
		for (int i = 0; i<m_rule.m_nth; ++i) {
            theta = i*dth;
            ksia = 0;
            ksib = sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-1)/sqrt(lC[0]*SQR(cos(theta))+lC[1]*SQR(sin(theta))-lC[2]);
            dksi = (ksib - ksia)/2;
            sksi = (ksib + ksia)/2;
			for (int j = 0; j<m_rule.m_nph; ++j) {
				ksi = sksi + dksi*m_rule.m_gp[j];
				wn = m_rule.m_gw[j] * dth*dksi;
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
double FEFiberIntegrationGauss::IntegratedFiberDensity()
{
    // establish local basis
    vec3d a[3], n0a;
    a[0] = vec3d(1,0,0);
    a[1] = vec3d(0,1,0);
    a[2] = vec3d(0,0,1);
    
    // initialize integrated fiber density
    double IFD = 1;
    double C = 0;
    
    double phi, theta;
    double ksia, ksib, dksi, sksi, ksi;
    double wn;
    double pi = 4*atan(1.0);
	double dth = 2 * pi / m_rule.m_nth;
    
    ksia = 0;
    ksib = 1;
    dksi = (ksib - ksia)/2;
    sksi = (ksib + ksia)/2;
    
    // loop over all integration points
	for (int i = 0; i<m_rule.m_nth; ++i) {
        theta = i*dth;
		for (int j = 0; j<m_rule.m_nph; ++j) {
			ksi = sksi + dksi*m_rule.m_gp[j];
			wn = m_rule.m_gw[j] * dth*dksi;
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
    return IFD;
}
*/
