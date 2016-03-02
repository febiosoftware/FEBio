//
//  FEFiberIntegrationGaussUC.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEFiberIntegrationGaussUC__
#define __FEBioMech__FEFiberIntegrationGaussUC__

#include "FEFiberIntegrationSchemeUC.h"

//----------------------------------------------------------------------------------
// Gauss integration scheme for continuous fiber distributions
//
class FEFiberIntegrationGaussUC : public FEFiberIntegrationSchemeUC
{
public:
    FEFiberIntegrationGaussUC(FEModel* pfem) : FEFiberIntegrationSchemeUC(pfem) { m_nph = 5; m_nth = 2*m_nph; }
    ~FEFiberIntegrationGaussUC() {}
	
	//! Initialization
	bool Init();
    
	//! Cauchy stress
	mat3ds DevStress(FEMaterialPoint& mp);
    
	// Spatial tangent
	tens4ds DevTangent(FEMaterialPoint& mp);
    
	//! Strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp);
    
    // Fiber density
    double IntegratedFiberDensity();
    
public:
	int             m_nph;	// number of gauss integration points along phi
    int             m_nth;  // number of trapezoidal integration points along theta
    vector<double>  m_gp;   // gauss points
    vector<double>  m_gw;   // gauss weights
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FEFiberIntegrationGaussUC__) */
