// FERigid.cpp: implementation of the FERigid class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigid.h"

// register the material with the framework
REGISTER_MATERIAL(FERigidMaterial, "rigid body");

// define the material parameters
BEGIN_PARAMETER_LIST(FERigidMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_PARAM_DOUBLE, "E");
	ADD_PARAMETER(m_v, FE_PARAM_DOUBLE, "v");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FERigid
//////////////////////////////////////////////////////////////////////

void FERigidMaterial::Init()
{
	FEElasticMaterial::Init();

	if (m_E <= 0) throw MaterialError("Invalid value for E");
	if (!INRANGE(m_v, -1.0, 0.5)) throw MaterialError("Invalid value for v");
}
