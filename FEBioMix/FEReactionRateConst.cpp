/*
 *  FEReactionRateConst.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 3/6/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

#include "FEReactionRateConst.h"

// register the material with the framework
REGISTER_MATERIAL(FEReactionRateConst, "constant");

// Material parameters for the FEMultiphasic material
BEGIN_PARAMETER_LIST(FEReactionRateConst, FEMaterial)
	ADD_PARAMETER(m_k, FE_PARAM_DOUBLE, "k");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEReactionRateConst::Init()
{
	FEMaterial::Init();

	if (m_k < 0) throw MaterialError("reaction rate constant must be positive");
	
}

