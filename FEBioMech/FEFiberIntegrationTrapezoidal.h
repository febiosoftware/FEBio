//
//  FEFiberIntegrationTrapezoidal.h
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 11/30/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#pragma once
#include "FEFiberIntegrationScheme.h"

//----------------------------------------------------------------------------------
// Trapezoidal integration scheme for 2D continuous fiber distributions
//
class FEFiberIntegrationTrapezoidal : public FEFiberIntegrationScheme
{
	class Iterator;

public:
    FEFiberIntegrationTrapezoidal(FEModel* pfem);
    ~FEFiberIntegrationTrapezoidal();

	// get iterator	
	FEFiberIntegrationSchemeIterator* GetIterator(FEMaterialPoint* mp) override;
    
private:
    int             m_nth;  // number of trapezoidal integration points along theta

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
