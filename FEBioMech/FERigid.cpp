// FERigid.cpp: implementation of the FERigid class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigid.h"
#include "FECore/FEModel.h"
#include "FECore/FERigidBody.h"

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

	m_binit = false;
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

	if (m_binit == false)
	{
		// get this rigid body's ID
		FERigidBody& rb = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(GetRigidBodyID()));

		// only set the rigid body com if this is the main rigid body material
		if (rb.GetMaterialID() == GetID()-1)
		{
			if (m_com == 1)
			{
				rb.SetCOM(m_rc);
			}
			else
			{
				rb.UpdateCOM();
			}
		}

		if (m_pmid  > -1)
		{
			FERigidMaterial* ppm = dynamic_cast<FERigidMaterial*>(GetFEModel()->GetMaterial(m_pmid-1));
			if (ppm == 0) throw MaterialError("parent of rigid material %s is not a rigid material\n", GetName());

			FERigidBody& prb = dynamic_cast<FERigidBody&>(*GetFEModel()->Object(ppm->GetRigidBodyID()));
			rb.m_prb = &prb;

			// we also need to open up all the RB's degree of freedoms
			m_bc[0] = 1;
			m_bc[1] = 1;
			m_bc[2] = 1;
			m_bc[3] = 1;
			m_bc[4] = 1;
			m_bc[5] = 1;
		}

		m_binit = true;
	}
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
