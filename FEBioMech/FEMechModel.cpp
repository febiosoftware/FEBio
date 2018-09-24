#include "stdafx.h"
#include "FEMechModel.h"
#include "FERigidSystem.h"
#include <FECore/FEMaterial.h>
#include "FERigidBody.h"
#include "FESolidMaterial.h"
#include "FERigidMaterial.h"

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
//! Initialize shells
void FEMechModel::InitShells()
{
	// Base class does most of the work
	FEModel::InitShells();

	// NOTE: This was moved here because I wanted to FEMaterial::IsRigid to FESolidMaterial::IsRigid
	//       This was part of the move to rid the FECore library of rigid stuff

	// Find the nodes that are on a non-rigid shell. 
	// These nodes will be assigned rotational degrees of freedom
	// TODO: Perhaps I should let the domains do this instead
	FEMesh& mesh = GetMesh();
	for (int i = 0; i<mesh.Nodes(); ++i) mesh.Node(i).m_nstate &= ~FENode::SHELL;
	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		if (dom.Class() == FE_DOMAIN_SHELL)
		{
			FESolidMaterial* pmat = dynamic_cast<FESolidMaterial*>(dom.GetMaterial());
			if (pmat && pmat->IsRigid() == false)
			{
				int N = dom.Elements();
				for (int i = 0; i<N; ++i)
				{
					FEElement& el = dom.ElementRef(i);
					int n = el.Nodes();
					for (int j = 0; j<n; ++j) mesh.Node(el.m_node[j]).m_nstate |= FENode::SHELL;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
// find a parameter value
FEParamValue FEMechModel::GetParameterValue(const ParamString& paramString)
{
	FEParamValue val = FEModel::GetParameterValue(paramString);

	if (val.isValid() == false)
	{
		// TODO: Move this to the rigid system
		ParamString next = paramString.next();
		if (next == "rigidbody")
		{
			FEMaterial* mat = 0;
			if (next.IDString()) mat = FindMaterial(next.IDString());
			if ((mat != 0) && (dynamic_cast<FERigidMaterial*>(mat)))
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

	return true;
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
