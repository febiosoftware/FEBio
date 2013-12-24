#pragma once
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
class FEMixDomainFactory : public FEDomainFactory
{
public:
	virtual int GetDomainType(const FE_Element_Spec& spec, FEMaterial* pmat);
	virtual FEDomain* CreateDomain(int dtype, FEMesh* pm, FEMaterial* pmat);
};
