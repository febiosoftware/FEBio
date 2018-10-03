//
//  FEFluidFSI.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 8/13/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FEFluidFSI.h"

#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidFSI, FEMaterial)
	// material properties
	ADD_PROPERTY(m_pSolid, "solid");
	ADD_PROPERTY(m_pFluid, "fluid");
END_FECORE_CLASS();

//============================================================================
// FEFSIMaterialPoint
//============================================================================
FEFSIMaterialPoint::FEFSIMaterialPoint(FEMaterialPoint* pt) : FEMaterialPoint(pt) {}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEFSIMaterialPoint::Copy()
{
    FEFSIMaterialPoint* pt = new FEFSIMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEFSIMaterialPoint::Serialize(DumpStream& ar)
{
    if (ar.IsSaving())
    {
        ar << m_w << m_aw << m_Jdot;
    }
    else
    {
        ar >> m_w >> m_aw >> m_Jdot;
    }
    
    FEMaterialPoint::Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEFSIMaterialPoint::Init()
{
    m_w = m_aw = vec3d(0,0,0);
    m_Jdot = 0;
    
    FEMaterialPoint::Init();
}

//============================================================================
// FEFluidFSI
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluidFSI constructor

FEFluidFSI::FEFluidFSI(FEModel* pfem) : FEMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEFluidFSI::CreateMaterialPointData()
{
    FEFluidMaterialPoint* fpt = new FEFluidMaterialPoint(m_pSolid->CreateMaterialPointData());
    return new FEFSIMaterialPoint(fpt);
}

//-----------------------------------------------------------------------------
// initialize
bool FEFluidFSI::Init()
{
    // set the solid density to zero (required for the solid of a FSI domain)
    m_pSolid->SetDensity(0.0);
    
    return FEMaterial::Init();
}
