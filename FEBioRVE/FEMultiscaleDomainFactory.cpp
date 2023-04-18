#include "stdafx.h"
#include "FEMultiscaleDomainFactory.h"
#include "FEMicroMaterial.h"
#include "FEMicroMaterial2O.h"
#include "FEElasticMultiscaleDomain1O.h"
#include "FEElasticMultiscaleDomain2O.h"

//==========================================================================================
FEDomain* FEMultiScaleDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel* pfem = pmat->GetFEModel();
	const char* sztype = 0;
	if      (dynamic_cast<FEMicroMaterial*    >(pmat)) sztype = "elastic-mm-solid";
	else if (dynamic_cast<FEMicroMaterial2O*  >(pmat)) sztype = "elastic-mm-solid2O";
	else if (dynamic_cast<FEElasticMaterial2O*>(pmat)) sztype = "elastic-solid2O";

	if (sztype)
	{
		FESolidDomain* pd = fecore_new<FESolidDomain>(sztype, pfem);
		if (pd) pd->SetMaterial(pmat);
		return pd;
	}
	else return 0;
}
