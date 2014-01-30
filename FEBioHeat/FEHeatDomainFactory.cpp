#include "FEHeatDomainFactory.h"
#include "FEHeatTransferMaterial.h"
#include "FEHeatSolidDomain.h"

//-----------------------------------------------------------------------------
int FEHeatDomainFactory::GetDomainType(const FE_Element_Spec& spec, FEMaterial* pmat)
{
	FE_Element_Shape eshape = spec.eshape;
	if (dynamic_cast<FEHeatTransferMaterial*>(pmat))
	{
		if ((eshape == ET_HEX8) || (eshape == ET_HEX20) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_TET15)) return FE_HEAT_SOLID_DOMAIN;
		else return 0;
	}
	return 0;
}

//-----------------------------------------------------------------------------
FEDomain* FEHeatDomainFactory::CreateDomain(int dtype, FEMesh* pm, FEMaterial* pmat)
{
	if (dtype == FE_HEAT_SOLID_DOMAIN) return new FEHeatSolidDomain(pm, pmat);
	return 0;
}
