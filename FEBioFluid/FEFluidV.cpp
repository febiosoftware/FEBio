//
//  FEFluidV.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/1/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#include "FEFluidV.h"
#include "FECore/FECoreKernel.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEFluidV, FEMaterial)
    // material properties
    ADD_PROPERTY(m_pFluid, "fluid");
END_FECORE_CLASS();

//============================================================================
// FEFluidVMaterialPoint
//============================================================================
FEFluidVMaterialPoint::FEFluidVMaterialPoint(FEMaterialPoint* pt) : FEMaterialPoint(pt) {}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEFluidVMaterialPoint::Copy()
{
    FEFluidVMaterialPoint* pt = new FEFluidVMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEFluidVMaterialPoint::Serialize(DumpStream& ar)
{
    if (ar.IsSaving())
    {
        ar << m_Jft << m_Jfp << m_dJft << m_dJfp;
    }
    else
    {
        ar >> m_Jft >> m_Jfp >> m_dJft >> m_dJfp;
    }
    
    FEMaterialPoint::Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEFluidVMaterialPoint::Init()
{
    m_Jft = m_Jfp = 1;
    m_dJft = m_dJfp = 0;
    
    FEMaterialPoint::Init();
}

//============================================================================
// FEFluidV
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluidV constructor

FEFluidV::FEFluidV(FEModel* pfem) : FEMaterial(pfem)
{
    m_pFluid = 0;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEFluidV::CreateMaterialPointData()
{
    FEFluidMaterialPoint* fpt = new FEFluidMaterialPoint();
    return new FEFluidVMaterialPoint(fpt);
}

//-----------------------------------------------------------------------------
// initialize
bool FEFluidV::Init()
{
    return FEMaterial::Init();
}
