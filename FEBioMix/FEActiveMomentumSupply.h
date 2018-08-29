//
//  FEActiveMomentumSupply.h
//  FEBioMix
//
//  Created by Gerard Ateshian on 12/9/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMix__FEActiveMomentumSupply__
#define __FEBioMix__FEActiveMomentumSupply__
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for active momentum supply.
//! These materials need to define the momentum supply and its tangents.
//!
class FEActiveMomentumSupply : public FEMaterial
{
public:
    FEActiveMomentumSupply(FEModel* pfem) : FEMaterial(pfem) {}
    virtual ~FEActiveMomentumSupply(){}
    
    //! active momentum supply
    virtual vec3d ActiveSupply(FEMaterialPoint& pt) = 0;
    
    //! tangent of active momentum supply with respect to strain
    virtual vec3d Tangent_ActiveSupply_Strain(FEMaterialPoint& mp) = 0;
    
    //! tangent of hydraulic permeability with respect to concentration
    vec3d Tangent_ActiveSupply_Concentration(FEMaterialPoint& mp, const int isol);
};


#endif /* defined(__FEBioMix__FEActiveMomentumSupply__) */
