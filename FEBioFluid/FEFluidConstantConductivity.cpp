//
//  FEFluidConstantConductivity.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/28/20.
//  Copyright © 2020 febio.org. All rights reserved.
//

#include "FEFluidConstantConductivity.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidConstantConductivity, FEFluidThermalConductivity)

    // material parameters
    ADD_PARAMETER(m_K   , FE_RANGE_GREATER_OR_EQUAL(0.0), "K");

END_FECORE_CLASS();

