/*
 *  FEReactionRateConst.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 3/6/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

#include "FEReactionRateConst.h"

// Material parameters for the FEReactionRateConst material
BEGIN_PARAMETER_LIST(FEReactionRateConst, FEMaterial)
	ADD_PARAMETER2(m_k, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "k");
END_PARAMETER_LIST();
