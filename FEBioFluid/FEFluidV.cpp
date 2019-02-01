//
//  FEFluidV.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/1/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#include "FEFluidV.h"
#include "FECore/FECoreKernel.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEFluidV, FEMaterial)
    // material properties
    ADD_PROPERTY(m_pFluid, "fluid");
END_FECORE_CLASS();

//============================================================================
// FEFluidV
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluidV constructor

FEFluidV::FEFluidV(FEModel* pfem) : FEMaterial(pfem)
{
    m_pFluid = 0;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEFluidV::CreateMaterialPointData()
{
    return new FEFluidMaterialPoint();
}

//-----------------------------------------------------------------------------
// initialize
bool FEFluidV::Init()
{
    return FEMaterial::Init();
}
