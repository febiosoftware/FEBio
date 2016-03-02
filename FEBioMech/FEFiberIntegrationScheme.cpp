//
//  FEFiberIntegrationScheme.cpp
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 11/16/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#include "FEFiberIntegrationScheme.h"

bool FEFiberIntegrationScheme::Init()
{
    if (m_pFmat->Init() == false) return false;
    if (m_pFDD->Init() == false) return false;
    
    // evaluate the integrated fiber density distribution
    if (m_pFDD->m_IFD == 1) m_pFDD->m_IFD = IntegratedFiberDensity();

	return true;
}
