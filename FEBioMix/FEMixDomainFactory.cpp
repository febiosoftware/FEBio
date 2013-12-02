#include "stdafx.h"
#include "FEMixDomainFactory.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include "FEBiphasicSolidDomain.h"
#include "FEBiphasicSoluteDomain.h"
#include "FETriphasicDomain.h"
#include "FEMultiphasicDomain.h"

//-----------------------------------------------------------------------------
int FEMixDomainFactory::GetDomainType(const FE_Element_Spec& spec, FEMaterial* pmat)
{
	FE_Element_Shape eshape = spec.eshape;
	if (dynamic_cast<FEBiphasic*>(pmat))
	{
		// biphasic elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_HEX20)) return FE_BIPHASIC_DOMAIN;
		else return 0;
	}
	if (dynamic_cast<FEBiphasicSolute*>(pmat))
	{
		// biphasic elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_HEX20)) return FE_BIPHASIC_SOLUTE_DOMAIN;
		else return 0;
	}
	else if (dynamic_cast<FETriphasic*>(pmat))
	{
		// triphasic elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_HEX20)) return FE_TRIPHASIC_DOMAIN;
		else return 0;
	}
	if (dynamic_cast<FEMultiphasic*>(pmat))
	{
		// multiphasic elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_HEX20)) return FE_MULTIPHASIC_DOMAIN;
		else return 0;
	}
	return 0;
}

//-----------------------------------------------------------------------------
FEDomain* FEMixDomainFactory::CreateDomain(int dtype, FEMesh* pm, FEMaterial* pmat)
{
	if (dtype == FE_BIPHASIC_DOMAIN       ) return new FEBiphasicSolidDomain (pm, pmat);
	if (dtype == FE_BIPHASIC_SOLUTE_DOMAIN) return new FEBiphasicSoluteDomain(pm, pmat);
	if (dtype == FE_TRIPHASIC_DOMAIN      ) return new FETriphasicDomain     (pm, pmat);
	if (dtype == FE_MULTIPHASIC_DOMAIN    ) return new FEMultiphasicDomain   (pm, pmat);
	return 0;
}
