//
//  FEFiberIntegrationGaussKronrodUC.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEFiberIntegrationGaussKronrodUC__
#define __FEBioMech__FEFiberIntegrationGaussKronrodUC__

#include "FEFiberIntegrationSchemeUC.h"

//----------------------------------------------------------------------------------
// Gauss integration scheme for continuous fiber distributions
//
class FEFiberIntegrationGaussKronrodUC : public FEFiberIntegrationSchemeUC
{
public:
    FEFiberIntegrationGaussKronrodUC(FEModel* pfem) : FEFiberIntegrationSchemeUC(pfem) { m_nph = 7; m_nth = 3; }
    ~FEFiberIntegrationGaussKronrodUC() {}
	
	//! Initialization
	void Init();
    
	//! Cauchy stress
	mat3ds DevStress(FEMaterialPoint& mp);
    
	// Spatial tangent
	tens4ds DevTangent(FEMaterialPoint& mp);
    
    // Fiber density
    void IntegratedFiberDensity(double& IFD);
    
public:
	int             m_nph;	// number of gauss integration points along phi
    int             m_nth;  // number of trapezoidal integration points along theta
    vector<double>  m_gp;   // gauss points
    vector<double>  m_gw;   // gauss weights
    bool            m_bfirst;
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FEFiberIntegrationGaussKronrodUC__) */
