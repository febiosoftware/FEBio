//
//  FEFiberIntegrationSchemeUC.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEFiberIntegrationSchemeUC__
#define __FEBioMech__FEFiberIntegrationSchemeUC__

#include "FEUncoupledMaterial.h"
#include "FEElasticFiberMaterialUC.h"
#include "FEFiberDensityDistribution.h"

//----------------------------------------------------------------------------------
// Base clase for integration schemes for continuous fiber distributions
//
class FEFiberIntegrationSchemeUC : public FEUncoupledMaterial
{
public:
    FEFiberIntegrationSchemeUC(FEModel* pfem) : FEUncoupledMaterial(pfem) {}
    
    bool Init();
    virtual void IntegratedFiberDensity(double& IFD) = 0;
    
public:
    FEElasticFiberMaterialUC*   m_pFmat;    // pointer to fiber material
    FEFiberDensityDistribution* m_pFDD;     // pointer to fiber density distribution
};


#endif /* defined(__FEBioMech__FEFiberIntegrationSchemeUC__) */
