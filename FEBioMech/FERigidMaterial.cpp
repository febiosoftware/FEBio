// FERigid.cpp: implementation of the FERigid class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigidMaterial.h"
#include "FECore/FEModel.h"
#include "FERigidSystem.h"
#include "FERigidBody.h"
#include "FEMechModel.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FERigidMaterial, FESolidMaterial)
	ADD_PARAMETER2(m_density, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "density"       );
	ADD_PARAMETER2(m_E      , FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "E"             );
	ADD_PARAMETER2(m_v      , FE_PARAM_DOUBLE, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v"             );
	ADD_PARAMETER(m_pmid   , FE_PARAM_INT   , "parent_id"     );
	ADD_PARAMETER(m_rc     , FE_PARAM_VEC3D , "center_of_mass");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// constructor
FERigidMaterial::FERigidMaterial(FEModel* pfem) : FESolidMaterial(pfem)
{
	m_com = 0;	// calculate COM automatically
	m_E = 1;
	m_v = 0;
	m_pmid = -1;

	m_binit = false;
}

//-----------------------------------------------------------------------------
void FERigidMaterial::SetParameter(FEParam& p)
{
	if (strcmp(p.name(), "center_of_mass") == 0)
	{
		m_com = 1;
	}
}

//-----------------------------------------------------------------------------
// Initialize rigid material data
bool FERigidMaterial::Init()
{
	if (FESolidMaterial::Init() == false) return false;

	if (m_binit == false)
	{
		// get this rigid body's ID
		FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
		FERigidSystem& rigid = *fem.GetRigidSystem();
		FERigidBody& rb = *rigid.Object(GetRigidBodyID());

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
			if (ppm == 0) return MaterialError("parent of rigid material %s is not a rigid material\n", GetName().c_str());

			FERigidBody& prb = *rigid.Object(ppm->GetRigidBodyID());
			rb.m_prb = &prb;

			// mark all degrees of freedom as prescribed
			rb.m_BC[0] = DOF_PRESCRIBED;
			rb.m_BC[1] = DOF_PRESCRIBED;
			rb.m_BC[2] = DOF_PRESCRIBED;
			rb.m_BC[3] = DOF_PRESCRIBED;
			rb.m_BC[4] = DOF_PRESCRIBED;
			rb.m_BC[5] = DOF_PRESCRIBED;
		}

		m_binit = true;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Serialize data to or from the dump file
void FERigidMaterial::Serialize(DumpStream &ar)
{
	// serialize base class parameters
	FESolidMaterial::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_com;
	}
	else
	{
		ar >> m_com;
	}
}
