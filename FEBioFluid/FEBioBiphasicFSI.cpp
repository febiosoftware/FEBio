//
//  FEBioBiphasicFSI.cpp
//  FEBioFluid
//
//  Created by Jay Shim on 12/30/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEBioBiphasicFSI.h"
#include <FECore/FECoreKernel.h>
#include "FEBiphasicFSISolver.h"
#include "FEBiphasicFSI.h"
#include "FEBiphasicFSIDomain3D.h"
#include "FEBiphasicFSIDomainFactory.h"
#include "FEBiphasicFSITraction.h"

//-----------------------------------------------------------------------------
const char* FEBioBiphasicFSI::GetVariableName(FEBioBiphasicFSI::BIPHASIC_FSI_VARIABLE var)
{
    switch (var)
    {
        case DISPLACEMENT                : return "displacement"               ; break;
        case VELOCITY                    : return "velocity"                   ; break;
        case SHELL_ROTATION              : return "shell rotation"             ; break;
        case SHELL_DISPLACEMENT          : return "shell displacement"         ; break;
        case SHELL_VELOCITY              : return "shell velocity"             ; break;
        case SHELL_ACCELERATION          : return "shell acceleration"         ; break;
        case RIGID_ROTATION              : return "rigid rotation"             ; break;
        case RELATIVE_FLUID_VELOCITY     : return "relative fluid velocity"    ; break;
        case RELATIVE_FLUID_ACCELERATION : return "relative fluid acceleration"; break;
        case FLUID_VELOCITY              : return "fluid velocity"             ; break;
        case FLUID_ACCELERATION          : return "fluid acceleration"         ; break;
        case FLUID_DILATATION            : return "fluid dilation"             ; break;
        case FLUID_DILATATION_TDERIV     : return "fluid dilation tderiv"      ; break;
    }
    assert(false);
    return nullptr;
}

void FEBioBiphasicFSI::InitModule()
{
    FECoreKernel& febio = FECoreKernel::GetInstance();
    
    // register domain
    febio.RegisterDomain(new FEBiphasicFSIDomainFactory);
    
    // define the fsi module
    febio.CreateModule("biphasic-FSI");
    febio.SetModuleDependency("fluid");
    
    REGISTER_FECORE_CLASS(FEBiphasicFSISolver, "biphasic-FSI");
    
    REGISTER_FECORE_CLASS(FEBiphasicFSI, "biphasic-FSI");
    
    REGISTER_FECORE_CLASS(FEBiphasicFSIDomain3D, "biphasic-FSI-3D");
    
    REGISTER_FECORE_CLASS(FEBiphasicFSITraction, "biphasic-FSI traction");
    
    febio.SetActiveModule(0);
}
