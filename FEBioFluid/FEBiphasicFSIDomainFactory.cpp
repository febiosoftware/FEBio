//
//  FEBiphasicFSIDomainFactory.cpp
//  FEBioFluid
//
//  Created by Jay Shim on 12/30/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEBiphasicFSIDomainFactory.h"
#include "FEBiphasicFSI.h"
#include <FECore/FEDomain.h>

//-----------------------------------------------------------------------------
FEDomain* FEBiphasicFSIDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
    FEModel* pfem = pmat->GetFEModel();
    FE_Element_Class eclass = spec.eclass;
    FE_Element_Shape eshape = spec.eshape;
    const char* sztype = 0;
    if (dynamic_cast<FEBiphasicFSI*>(pmat))
    {
        // fluid elements
        if (eclass == FE_ELEM_SOLID) sztype = "biphasic-FSI-3D";
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
