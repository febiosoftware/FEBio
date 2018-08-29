//
//  FEActiveMomentumSupply.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 12/9/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEActiveMomentumSupply.h"

//-----------------------------------------------------------------------------
// Derivative of active momentum supply w.r.t. solute concentration at material point
// Set this to zero by default because biphasic problems do not require it
vec3d FEActiveMomentumSupply::Tangent_ActiveSupply_Concentration(FEMaterialPoint& pt, const int isol)
{
    return vec3d(0,0,0);
}

