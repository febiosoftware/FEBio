//
//  FEBioBiphasicFSI.hpp
//  FEBioFluid
//
//  Created by Jay Shim on 12/30/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#pragma once
#include "febiofluid_api.h"

namespace FEBioBiphasicFSI {
    
    FEBIOFLUID_API void InitModule();
    
    enum BIPHASIC_FSI_VARIABLE {
        DISPLACEMENT,
        VELOCITY,
        SHELL_ROTATION,
        SHELL_DISPLACEMENT,
        SHELL_VELOCITY,
        SHELL_ACCELERATION,
        RIGID_ROTATION,
        RELATIVE_FLUID_VELOCITY,
        RELATIVE_FLUID_ACCELERATION,
        FLUID_VELOCITY,
        FLUID_ACCELERATION,
        FLUID_DILATATION,
        FLUID_DILATATION_TDERIV
    };
    
    FEBIOFLUID_API const char* GetVariableName(BIPHASIC_FSI_VARIABLE var);
}
