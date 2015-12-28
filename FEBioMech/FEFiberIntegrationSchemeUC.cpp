//
//  FEFiberIntegrationSchemeUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEFiberIntegrationSchemeUC.h"

bool FEFiberIntegrationSchemeUC::Init()
{
    if (m_pFmat->Init() == false) return false;
    if (m_pFDD->Init() == false) return false;
    
    // evaluate the integrated fiber density distribution
    IntegratedFiberDensity(m_pFDD->m_IFD);

	return true;
}
