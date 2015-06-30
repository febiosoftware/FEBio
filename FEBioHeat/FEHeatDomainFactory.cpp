#include "FEHeatDomainFactory.h"
#include "FEHeatTransferMaterial.h"
#include "FEHeatSolidDomain.h"
#include "FEThermoElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! Use the material and the element type to determine the domain type.
int FEHeatDomainFactory::GetDomainType(const FE_Element_Spec& spec, FEMaterial* pmat)
{
	FE_Element_Shape eshape = spec.eshape;
	if (dynamic_cast<FEHeatTransferMaterial*>(pmat))
	{
		if ((eshape == ET_HEX8) || (eshape == ET_HEX20) || (eshape == ET_HEX27) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_TET15)) return FE_HEAT_SOLID_DOMAIN;
		else return 0;
	}
	if (dynamic_cast<FEThermoElasticMaterial*>(pmat))
	{
		if ((eshape == ET_HEX8) || (eshape == ET_HEX20) || (eshape == ET_HEX27) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_TET15)) return FE_THERMOELASTIC_DOMAIN;
		else return 0;
	}
	return 0;
}

//-----------------------------------------------------------------------------
//! Create a domain from the given domain type (dtype)
FEDomain* FEHeatDomainFactory::CreateDomain(int dtype, FEMesh* pm, FEMaterial* pmat)
{
	if (dtype == FE_HEAT_SOLID_DOMAIN   ) return new FEHeatSolidDomain(pm, pmat);
	if (dtype == FE_THERMOELASTIC_DOMAIN) return new FEThermoElasticSolidDomain(pm, pmat);
	return 0;
}
