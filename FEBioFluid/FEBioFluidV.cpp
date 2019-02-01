//
//  FEBioFluidV.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 1/30/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#include "FEBioFluidV.h"
#include <FECore/FECoreKernel.h>
#include "FEFluidVSolver.h"
#include "FEFluidVDomain3D.h"
#include "FEFluidVDomainFactory.h"

void FEBioFluidV::InitModule()
{
    FECoreKernel& febio = FECoreKernel::GetInstance();
    
    // register domain
    febio.RegisterDomain(new FEFluidVDomainFactory);
    
    // define the fluidV module
    febio.CreateModule("fluidV");
    febio.SetModuleDependency("fluid");
    
    REGISTER_FECORE_CLASS(FEFluidVSolver, "fluidV");
    
    REGISTER_FECORE_CLASS(FEFluidVDomain3D, "fluidV-3D");
    
    febio.SetActiveModule(0);
}
