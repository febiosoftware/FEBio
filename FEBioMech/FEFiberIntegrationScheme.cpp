//
//  FEFiberIntegrationScheme.cpp
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 11/16/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#include "FEFiberIntegrationScheme.h"

void FEFiberIntegrationScheme::Init()
{
    m_pFmat->Init();
    m_pFDD->Init();
    
    // evaluate the integrated fiber density distribution
    IntegratedFiberDensity(m_pFDD->m_IFD);
}
