//
//  FEMembraneReactionRateConst.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 3/6/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#include "FEMembraneReactionRateConst.h"

// Material parameters for the FEMembraneReactionRateConst material
BEGIN_PARAMETER_LIST(FEMembraneReactionRateConst, FEMembraneReactionRate)
	ADD_PARAMETER(m_k, FE_RANGE_GREATER_OR_EQUAL(0.0), "k");
END_PARAMETER_LIST();
