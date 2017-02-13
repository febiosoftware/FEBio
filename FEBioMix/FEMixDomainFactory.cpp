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
FEDomain* FEMixDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel* pfem = pmat->GetFEModel();
	FE_Element_Class eclass = spec.eclass;
	const char* sztype = 0;
	if (dynamic_cast<FEBiphasic*>(pmat))
	{
		// biphasic elements
		if (eclass == FE_ELEM_SOLID) sztype = "biphasic-solid";
        else if (eclass == FE_ELEM_SHELL) sztype = "biphasic-shell";
		else return 0;
	}
	if (dynamic_cast<FEBiphasicSolute*>(pmat))
	{
		// biphasic solute elements
		if (eclass == FE_ELEM_SOLID) sztype = "biphasic-solute-solid";
        else if (eclass == FE_ELEM_SHELL) sztype = "biphasic-solute-shell";
		else return 0;
	}
	else if (dynamic_cast<FETriphasic*>(pmat))
	{
		// triphasic elements
		if (eclass == FE_ELEM_SOLID) sztype = "triphasic-solid";
		else return 0;
	}
	if (dynamic_cast<FEMultiphasic*>(pmat))
	{
		// multiphasic elements
		if (eclass == FE_ELEM_SOLID)  sztype = "multiphasic-solid";
        else if (eclass == FE_ELEM_SHELL) sztype = "multiphasic-shell";
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
