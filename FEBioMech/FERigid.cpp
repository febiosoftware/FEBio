// FERigid.cpp: implementation of the FERigid class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigid.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FERigidMaterial, FESolidMaterial)
	ADD_PARAMETER(m_density, FE_PARAM_DOUBLE, "density"       );
	ADD_PARAMETER(m_E      , FE_PARAM_DOUBLE, "E"             );
	ADD_PARAMETER(m_v      , FE_PARAM_DOUBLE, "v"             );
	ADD_PARAMETER(m_pmid   , FE_PARAM_INT   , "parent_id"     );
	ADD_PARAMETER(m_rc     , FE_PARAM_VEC3D , "center_of_mass");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// constructor
FERigidMaterial::FERigidMaterial(FEModel* pfem) : FESolidMaterial(pfem)
{
	m_com = 0;	// calculate COM automatically
	for (int i=0; i<6; ++i)
	{
		m_bc[i] =  0;	// rigid bodies are initially free
		m_fc[i] = -1;
		m_fs[i] =  0;
	}
	m_E = 1;
	m_v = 0;
	m_pmid = -1;
}

//-----------------------------------------------------------------------------
void FERigidMaterial::SetParameter(FEParam& p)
{
	if (strcmp(p.m_szname, "center_of_mass") == 0)
	{
		m_com = 1;
	}
}

//-----------------------------------------------------------------------------
// Initialize rigid material data
void FERigidMaterial::Init()
{
	FESolidMaterial::Init();

	if (m_E <= 0) throw MaterialError("Invalid value for E");
	if (!IN_RIGHT_OPEN_RANGE(m_v, -1.0, 0.5)) throw MaterialError("Invalid value for v");
}

//-----------------------------------------------------------------------------
//! Serialize data to or from the dump file
void FERigidMaterial::Serialize(DumpFile &ar)
{
	// serialize base class parameters
	FESolidMaterial::Serialize(ar);

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
