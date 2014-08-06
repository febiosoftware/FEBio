//
//  FEFiberIntegrationSchemeUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEFiberIntegrationSchemeUC.h"

void FEFiberIntegrationSchemeUC::Init()
{
    m_pFmat->Init();
    m_pFDD->Init();
    
    // evaluate the integrated fiber density distribution
    IntegratedFiberDensity(m_pFDD->m_IFD);
}
