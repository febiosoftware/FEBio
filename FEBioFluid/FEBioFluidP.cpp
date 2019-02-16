//
//  FEBioFluidP.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/16/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#include "FEBioFluidP.h"
#include <FECore/FECoreKernel.h>
#include "FEFluidSolver.h"
#include "FEFluidPDomain3D.h"
#include "FEFluidPDomainFactory.h"
#include "FEFluidPResistanceBC.h"

void FEBioFluidP::InitModule()
{
    FECoreKernel& febio = FECoreKernel::GetInstance();
    
    // register domain
    febio.RegisterDomain(new FEFluidPDomainFactory);
    
    // define the fluidP module
    febio.CreateModule("fluidP");
    febio.SetModuleDependency("fluid");

    REGISTER_FECORE_CLASS(FEFluidPDomain3D, "fluidP-3D");
    
    REGISTER_FECORE_CLASS(FEFluidPResistanceBC, "fluidP resistance");
    
    febio.SetActiveModule(0);
}
