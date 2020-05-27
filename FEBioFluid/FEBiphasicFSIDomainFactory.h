//
//  FEBiphasicFSIDomainFactory.hpp
//  FEBioFluid
//
//  Created by Jay Shim on 12/30/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#pragma once
#include <FECore/FECoreKernel.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FEBiphasicFSIDomainFactory : public FEDomainFactory
{
public:
    virtual FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
