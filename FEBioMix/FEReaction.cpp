//
//  FEReaction.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 3/4/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEReaction.h"
#include "FECore/FEElementTraits.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include <FECore/fecore_error.h>
#include "FEMultiphasic.h"
#include <stdlib.h>

//-----------------------------------------------------------------------------
FEReaction::FEReaction(FEModel* pfem) : FEMaterial(pfem)
{
    m_pMP = 0;
}

//-----------------------------------------------------------------------------
bool FEReaction::Init()
{
    // make sure the parent class is set
    assert(m_pMP);
    if (m_pMP == 0) return fecore_error("Parent class not set");
    
    // now call base class
    if (FEMaterial::Init() == false) return false;
    
    return true;
}

