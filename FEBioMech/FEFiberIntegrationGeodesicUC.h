//
//  FEFiberIntegrationGeodesicUC.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEFiberIntegrationGeodesicUC__
#define __FEBioMech__FEFiberIntegrationGeodesicUC__

#include "FEFiberIntegrationSchemeUC.h"
#include "geodesic.h"

//----------------------------------------------------------------------------------
// Geodesic dome integration scheme for continuous fiber distributions
//
class FEFiberIntegrationGeodesicUC : public FEFiberIntegrationSchemeUC
{
public:
    FEFiberIntegrationGeodesicUC(FEModel* pfem) : FEFiberIntegrationSchemeUC(pfem) { m_nres = 0; }
    ~FEFiberIntegrationGeodesicUC() {}
	
	//! Initialization
	void Init();
    
	//! Cauchy stress
	mat3ds DevStress(FEMaterialPoint& mp);
    
	// Spatial tangent
	tens4ds DevTangent(FEMaterialPoint& mp);
    
	//! Strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp);
    
    // Fiber density
    void IntegratedFiberDensity(double& IFD);
    
public:
	int             m_nres;	// resolution
    int             m_nint; // number of integration points
	double          m_cth[NSTH];
    double          m_sth[NSTH];
    double          m_cph[NSTH];
    double          m_sph[NSTH];
    double          m_w[NSTH];
    bool            m_bfirst;
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FEFiberIntegrationGeodesicUC__) */
