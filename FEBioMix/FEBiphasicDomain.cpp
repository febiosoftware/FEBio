//
//  FEBiphasicDomain.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 12/1/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "FEBiphasicDomain.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FEBiphasicDomain::FEBiphasicDomain(FEModel* pfem) : FEElasticDomain(pfem)
{
    m_pMat = 0;
    m_dofP = pfem->GetDOFIndex("p");
    m_dofQ = pfem->GetDOFIndex("q");
    
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
}
