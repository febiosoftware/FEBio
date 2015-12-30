//
//  FEFiberIntegrationTrapezoidalUC.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEFiberIntegrationTrapezoidalUC__
#define __FEBioMech__FEFiberIntegrationTrapezoidalUC__

#include "FEFiberIntegrationSchemeUC.h"

//----------------------------------------------------------------------------------
// Trapezoidal integration scheme for 2D continuous fiber distributions
//
class FEFiberIntegrationTrapezoidalUC : public FEFiberIntegrationSchemeUC
{
public:
    FEFiberIntegrationTrapezoidalUC(FEModel* pfem) : FEFiberIntegrationSchemeUC(pfem) { m_nth = 12; }
    ~FEFiberIntegrationTrapezoidalUC() {}
	
	//! Cauchy stress
	mat3ds DevStress(FEMaterialPoint& mp);
    
	// Spatial tangent
	tens4ds DevTangent(FEMaterialPoint& mp);
    
	//! Strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp);
    
    // Fiber density
    void IntegratedFiberDensity(double& IFD);
    
public:
    int             m_nth;  // number of trapezoidal integration points along theta
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FEFiberIntegrationTrapezoidalUC__) */
