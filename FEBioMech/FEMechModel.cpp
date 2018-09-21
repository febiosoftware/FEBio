#include "stdafx.h"
#include "FEMechModel.h"
#include "FERigidSystem.h"
#include <FECore/FEMaterial.h>
#include "FERigidBody.h"

//-----------------------------------------------------------------------------
FEMechModel::FEMechModel()
{
	// create a rigid system
	m_prs = new FERigidSystem(this);
}

//-----------------------------------------------------------------------------
// get the rigid system
FERigidSystem* FEMechModel::GetRigidSystem() { return m_prs; }

//-----------------------------------------------------------------------------
// clear all model data
void FEMechModel::Clear()
{
	m_prs->Clear();

	FEModel::Clear();
}

//-----------------------------------------------------------------------------
bool FEMechModel::InitRigidSystem()
{
	return m_prs->Init();
}

//-----------------------------------------------------------------------------
// model activation
void FEMechModel::Activate()
{
	FEModel::Activate();

	// activate rigid components
	m_prs->Activate();
}

//-----------------------------------------------------------------------------
bool FEMechModel::Reset()
{
	if (FEModel::Reset() == false) return false;
	return m_prs->Reset();
}

//-----------------------------------------------------------------------------
// Defined in FEModel.cpp
FEParamValue GetParameterComponent(const ParamString& paramName, FEParam* param);

//-----------------------------------------------------------------------------
// find a parameter value
FEParamValue FEMechModel::GetParameterValue(const ParamString& paramString)
{
	FEParamValue val = FEModel::GetParameterValue(paramString);

	if (val.isValid() == false)
	{
		ParamString next = paramString.next();
		if (next == "rigidbody")
		{
			FEMaterial* mat = 0;
			if (next.IDString()) mat = FindMaterial(next.IDString());
			if ((mat != 0) && (mat->IsRigid()))
			{
				ParamString paramName = next.next();

				// the rigid bodies are dealt with differently
				int nmat = mat->GetID() - 1;
				int NRB = m_prs->Objects();
				for (int i = 0; i<NRB; ++i)
				{
					FERigidBody* ob = m_prs->Object(i);
					if (ob && (ob->GetMaterialID() == nmat))
					{
						FEParam* pp = ob->FindParameter(paramName);
						return GetParameterComponent(paramName.last(), pp);
					}
				}
			}
		}
	}

	return val;
}

//-----------------------------------------------------------------------------
//! evaluate all parameter lists
bool FEMechModel::EvaluateAllParameterLists()
{
	if (FEModel::EvaluateAllParameterLists() == false) return false;

	// give the rigid system a chance
	if (m_prs->EvaluateParameterLists() == false) return false;
}

//-----------------------------------------------------------------------------
//! find a model componnet from its class ID
FEModelComponent* FEMechModel::FindModelComponent(int nid)
{
	FEModelComponent* mc = FEModel::FindModelComponent(nid);
	if (mc == 0) mc = m_prs->FindModelComponent(nid);
	return mc;
}

//-----------------------------------------------------------------------------
//! serialize data for restarts
void FEMechModel::Serialize(DumpStream& ar)
{
	FEModel::Serialize(ar);

	// stream rigid body data
	m_prs->Serialize(ar);
	ar.check();
}

//-----------------------------------------------------------------------------
//! Build the matrix profile for this model
void FEMechModel::BuildMatrixProfile(FEGlobalMatrix& G, bool breset)
{
	FEModel::BuildMatrixProfile(G, breset);

	if (breset)
	{
		// Add rigid bodies to the profile
		m_prs->BuildMatrixProfile(G);
	}
}
