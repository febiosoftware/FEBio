#pragma once
#include <FECore/FECoreKernel.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FEFluidPDomainFactory : public FEDomainFactory
{
public:
    virtual FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
