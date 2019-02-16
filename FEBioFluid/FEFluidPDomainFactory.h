//
//  FEFluidPDomainFactory.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/16/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#ifndef FEFluidPDomainFactory_hpp
#define FEFluidPDomainFactory_hpp

#include <FECore/FECoreKernel.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FEFluidPDomainFactory : public FEDomainFactory
{
public:
    virtual FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};

#endif /* FEFluidPDomainFactory_hpp */
