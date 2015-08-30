#pragma once
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
class FEMixDomainFactory : public FEDomainFactory
{
public:
	virtual FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
