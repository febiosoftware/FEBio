/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEMechModel.h"
#include "FERigidSystem.h"
#include <FECore/FEMaterial.h>
#include <FECore/FEDomain.h>
#include "FERigidBody.h"
#include "FESolidMaterial.h"
#include "FERigidMaterial.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEMechModel, FEModel)
	ADD_PROPERTY(m_prs->RigidBodyList(), "rigid_body");
END_FECORE_CLASS();


//-----------------------------------------------------------------------------
FEMechModel::FEMechModel()
{
	// create a rigid system
	m_prs = new FERigidSystem(this);
}

//-----------------------------------------------------------------------------
// get the rigid system
FERigidSystem* FEMechModel::GetRigidSystem() 
{ 
	return m_prs; 
}

//-----------------------------------------------------------------------------
// clear all model data
void FEMechModel::Clear()
{
	m_prs->Clear();

	FEModel::Clear();
}

//-----------------------------------------------------------------------------
// number of rigid bodies
int FEMechModel::RigidBodies() const
{
	return m_prs->Objects();
}

//-----------------------------------------------------------------------------
bool FEMechModel::InitRigidSystem()
{
	return m_prs->Init();
}

//-----------------------------------------------------------------------------
// get a rigid body
FERigidBody* FEMechModel::GetRigidBody(int n)
{
	return m_prs->Object(n);
}

//-----------------------------------------------------------------------------
// find a rigid body from a material ID
int FEMechModel::FindRigidbodyFromMaterialID(int matId)
{
	int NRB = RigidBodies();
	for (int i = 0; i<NRB; ++i)
	{
		FERigidBody& rb = *GetRigidBody(i);
		if (rb.GetMaterialID() == matId) return i;
	}
	return -1;
}

//-----------------------------------------------------------------------------
// return number or rigid prescribed BCs
int FEMechModel::RigidPrescribedBCs() const
{
	return m_prs->PrescribedBCs();
}

//-----------------------------------------------------------------------------
// return the rigid prescribed displacement
FERigidPrescribedBC* FEMechModel::GetRigidPrescribedBC(int i)
{
	return m_prs->PrescribedBC(i);
}

//-----------------------------------------------------------------------------
// add a rigid presribed BC
void FEMechModel::AddRigidPrescribedBC(FERigidPrescribedBC* pDC)
{
	m_prs->AddPrescribedBC(pDC);
}

//-----------------------------------------------------------------------------
// add a rigid fixed BC
void FEMechModel::AddRigidFixedBC(FERigidFixedBC* pBC)
{
	m_prs->AddFixedBC(pBC);
}

//-----------------------------------------------------------------------------
// add a rigid initial condition
void FEMechModel::AddRigidInitialCondition(FERigidIC* pIC)
{
	m_prs->AddInitialCondition(pIC);
}

//-----------------------------------------------------------------------------
// model activation
void FEMechModel::Activate()
{
    // activate rigid components
    m_prs->Activate();
    
	FEModel::Activate();
}

//-----------------------------------------------------------------------------
void FEMechModel::Reactivate()
{
	FEModel::Reactivate();
	m_prs->Activate();
}

//-----------------------------------------------------------------------------
bool FEMechModel::Reset()
{
	if (m_prs->Reset() == false) return false;
	return FEModel::Reset();
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
			if ((pmat == 0) || (pmat->IsRigid() == false))
			{
				int N = dom.Elements();
				for (int i = 0; i<N; ++i)
				{
					FEElement& el = dom.ElementRef(i);
					int n = el.Nodes();
					for (int j = 0; j < n; ++j)
					{
						mesh.Node(el.m_node[j]).m_nstate |= FENode::SHELL;

						// TODO: Not sure why, but it looks like this is needed otherwise
						//       the rigid bodies lock. 
						mesh.Node(el.m_node[j]).m_nstate |= FENode::RIGID_CLAMP;
					}
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
		return m_prs->GetParameterValue(paramString);
	}

	return val;
}

//-----------------------------------------------------------------------------
void FEMechModel::SerializeGeometry(DumpStream& ar)
{
	FEModel::SerializeGeometry(ar);
	m_prs->Serialize(ar);
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

//-----------------------------------------------------------------------------
// update rigid part of mesh
void FEMechModel::UpdateRigidMesh()
{
	m_prs->UpdateMesh();
}
