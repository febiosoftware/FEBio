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


//-----------------------------------------------------------------------------
// Initialize rigid material data
void FERigidMaterial::Init()
{
	FEElasticMaterial::Init();

	if (m_E <= 0) throw MaterialError("Invalid value for E");
	if (!IN_RIGHT_OPEN_RANGE(m_v, -1.0, 0.5)) throw MaterialError("Invalid value for v");
}

//-----------------------------------------------------------------------------
//! Serialize data to or from the dump file
void FERigidMaterial::Serialize(DumpFile &ar)
{
	// serialize base class parameters
	FEElasticMaterial::Serialize(ar);

	// TODO: do we really need to store this data?
	if (ar.IsSaving())
	{
		ar.write(m_bc, sizeof(int), 6);
		ar.write(m_fc, sizeof(int), 6);
		ar.write(m_fs, sizeof(double), 6);
	}
	else
	{
		ar.read(m_bc, sizeof(int), 6);
		ar.read(m_fc, sizeof(int), 6);
		ar.read(m_fs, sizeof(double), 6);
	}
}
