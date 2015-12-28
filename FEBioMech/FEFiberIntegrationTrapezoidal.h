//
//  FEFiberIntegrationTrapezoidal.h
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 11/30/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#pragma once
#include "FEFiberIntegrationScheme.h"

//----------------------------------------------------------------------------------
// Trapezoidal integration scheme for 2D continuous fiber distributions
//
class FEFiberIntegrationTrapezoidal : public FEFiberIntegrationScheme
{
public:
    FEFiberIntegrationTrapezoidal(FEModel* pfem) : FEFiberIntegrationScheme(pfem) { m_nth = 12; }
    ~FEFiberIntegrationTrapezoidal() {}
	
	//! Initialization
	bool Init();
    
	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp);
    
	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp);
    
	//! Strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp);
    
    // Fiber density
    void IntegratedFiberDensity(double& IFD);
    
public:
    int             m_nth;  // number of trapezoidal integration points along theta

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
