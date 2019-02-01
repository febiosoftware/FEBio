//
//  FEFluidVDomainFactory.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 1/30/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#ifndef FEFluidVDomainFactory_hpp
#define FEFluidVDomainFactory_hpp

#include <FECore/FECoreKernel.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FEFluidVDomainFactory : public FEDomainFactory
{
public:
    virtual FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};

#endif /* FEFluidVDomainFactory_hpp */
