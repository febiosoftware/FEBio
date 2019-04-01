/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
			if ((pmat == 0) || (pmat->IsRigid() == false))
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
		return m_prs->GetParameterValue(paramString);
	}

	return val;
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
