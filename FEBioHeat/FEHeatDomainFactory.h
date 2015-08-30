#pragma once
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
class FEHeatDomainFactory : public FEDomainFactory
{
public:
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
