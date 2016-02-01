#include "stdafx.h"
#include "FEFluidDomainFactory.h"
#include "FEFluid.h"
#include "FEFluidDomain.h"
#include <FECore/FEDomain.h>

//-----------------------------------------------------------------------------
FEDomain* FEFluidDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel* pfem = pmat->GetFEModel();
	FE_Element_Shape eshape = spec.eshape;
	const char* sztype = 0;
	if (dynamic_cast<FEFluid*>(pmat))
	{
		// fluid elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_HEX20) || (eshape == ET_HEX27)) sztype = "fluid";
        else if ((eshape == ET_QUAD4) || (eshape == ET_QUAD8) || (eshape == ET_QUAD9) || (eshape == ET_TRI3) || (eshape == ET_TRI6)) sztype = "fluid2D";
		else return 0;
	}

	if (sztype)
	{
		FEDomain* pd = fecore_new<FEDomain>(FEDOMAIN_ID, sztype, pfem);
		if (pd) pd->SetMaterial(pmat);
		return pd;
	}
	else return 0;
}
