#include "stdafx.h"
#include "FEFluidDomainFactory.h"
#include "FEFluid.h"
#include "FEFluidDomain.h"
#include <FECore/FEDomain.h>

//-----------------------------------------------------------------------------
FEDomain* FEFluidDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel* pfem = pmat->GetFEModel();
	FE_Element_Class eclass = spec.eclass;
	FE_Element_Shape eshape = spec.eshape;
	const char* sztype = 0;
	if (dynamic_cast<FEFluid*>(pmat))
	{
		// fluid elements
		if      (eclass==FE_ELEM_SOLID) sztype = "fluid";
        else if (eclass==FE_ELEM_2D   ) sztype = "fluid2D";
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
