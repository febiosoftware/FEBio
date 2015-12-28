//
//  FEFiberIntegrationScheme.h
//
//  Created by Gerard Ateshian on 11/16/13.
//

#pragma once
#include "FEElasticMaterial.h"
#include "FEElasticFiberMaterial.h"
#include "FEFiberDensityDistribution.h"

//----------------------------------------------------------------------------------
// Base clase for integration schemes for continuous fiber distributions
//
class FEFiberIntegrationScheme : public FEElasticMaterial
{
public:
    FEFiberIntegrationScheme(FEModel* pfem) : FEElasticMaterial(pfem) {}
    
    bool Init();
    virtual void IntegratedFiberDensity(double& IFD) = 0;

public:
    FEElasticFiberMaterial*     m_pFmat;    // pointer to fiber material
    FEFiberDensityDistribution* m_pFDD;     // pointer to fiber density distribution
};
