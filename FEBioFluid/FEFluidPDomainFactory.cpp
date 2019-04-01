#include "stdafx.h"
#include "FEFluidPDomainFactory.h"
#include "FEFluidP.h"
#include <FECore/FEDomain.h>

//-----------------------------------------------------------------------------
FEDomain* FEFluidPDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
    FEModel* pfem = pmat->GetFEModel();
    FE_Element_Class eclass = spec.eclass;
    FE_Element_Shape eshape = spec.eshape;
    const char* sztype = 0;
    if (dynamic_cast<FEFluidP*>(pmat))
    {
        // fluid elements
        if (eclass == FE_ELEM_SOLID) sztype = "fluidP-3D";
        else return 0;
    }
    
    if (sztype)
    {
        FEDomain* pd = fecore_new<FEDomain>(sztype, pfem);
        if (pd) pd->SetMaterial(pmat);
        return pd;
    }
    else return 0;
}
