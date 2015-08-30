#include "FEHeatDomainFactory.h"
#include "FEHeatTransferMaterial.h"
#include "FEHeatSolidDomain.h"
#include "FEThermoElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! Use the material and the element type to determine the domain type.
FEDomain* FEHeatDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel& fem = *pmat->GetFEModel();
	FE_Element_Shape eshape = spec.eshape;
	if (dynamic_cast<FEHeatTransferMaterial*>(pmat))
	{
		switch (eshape)
		{
		case ET_HEX8:
		case ET_HEX20:
		case ET_HEX27:
		case ET_PENTA6:
		case ET_TET4:
		case ET_TET10:
		case ET_TET15:
			{
				FEDomain* pd = new FEHeatSolidDomain(&fem);
				pd->SetMaterial(pmat);
				return pd;
			}
			break;
		}
	}
	if (dynamic_cast<FEThermoElasticMaterial*>(pmat))
	{
		switch (eshape)
		{
		case ET_HEX8:
		case ET_HEX20:
		case ET_HEX27:
		case ET_PENTA6:
		case ET_TET4:
		case ET_TET10:
		case ET_TET15:
			{
				FEDomain* pd = new FEThermoElasticSolidDomain(&fem);
				pd->SetMaterial(pmat);
				return pd;
			}
		}
	}
	return 0;
}
