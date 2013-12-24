#pragma once
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
class FEHeatDomainFactory : public FEDomainFactory
{
public:
	int GetDomainType(const FE_Element_Spec& spec, FEMaterial* pmat);
	FEDomain* CreateDomain(int dtype, FEMesh* pm, FEMaterial* pmat);
};
