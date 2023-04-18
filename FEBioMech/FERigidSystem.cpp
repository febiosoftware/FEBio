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
#include "FERigidSystem.h"
#include "FERigidBody.h"
#include <FECore/FEModel.h>
#include <FECore/FEMaterial.h>
#include <FECore/FEDomain.h>
#include "RigidBC.h"
#include <FECore/FEGlobalMatrix.h>
#include "FERigidMaterial.h"

//-----------------------------------------------------------------------------
//! constructor
FERigidSystem::FERigidSystem(FEModel* pfem) : m_fem(*pfem)
{
}

//-----------------------------------------------------------------------------
//! Add a rigid body
void FERigidSystem::AddRigidBody(FERigidBody* prb)
{
	if (prb) m_RB.push_back(prb);
}

//-----------------------------------------------------------------------------
FERigidSystem::~FERigidSystem()
{
	Clear();
}

//-----------------------------------------------------------------------------
int FERigidSystem::Objects() const
{
	return (int) m_RB.size();
}

//-----------------------------------------------------------------------------
std::vector<FERigidBody*>& FERigidSystem::RigidBodyList()
{
	return m_RB;
}

//-----------------------------------------------------------------------------
//! Get a rigid body
FERigidBody* FERigidSystem::Object(int i)
{
	return m_RB[i];
}

//-----------------------------------------------------------------------------
//! delete all rigid bodies
void FERigidSystem::Clear()
{
	for (int i=0; i<(int)m_RB.size (); ++i) delete m_RB [i]; m_RB.clear ();
	for (int i=0; i<(int)m_RBC.size(); ++i) delete m_RBC[i]; m_RBC.clear();
	for (int i=0; i<(int)m_RDC.size(); ++i) delete m_RDC[i]; m_RDC.clear();
	for (int i=0; i<(int)m_RIC.size(); ++i) delete m_RIC[i]; m_RIC.clear();
}

//-----------------------------------------------------------------------------
void FERigidSystem::Serialize(DumpStream& ar)
{
	if (ar.IsShallow())
	{
		for (int i=0; i<(int) m_RB.size(); ++i) m_RB[i]->Serialize(ar);
	}
	else
	{
		if (ar.IsSaving())
		{
			// rigid objects
			int nrb = Objects();
			ar << nrb;
			for (int i=0; i<nrb; ++i) Object(i)->Serialize(ar);

			// fixed rigid body dofs
			ar & m_RBC;

			// rigid body displacements
			ar & m_RDC;

			// rigid body initial conditions
			ar & m_RIC;
		}
		else
		{
			Clear();

			// rigid bodies
			int nrb = 0;
			ar >> nrb;
			for (int i=0; i<nrb; ++i)
			{
				FERigidBody* prb = new FERigidBody(&m_fem);
				prb->Serialize(ar);
				AddRigidBody(prb);
			}

			// fixed rigid body dofs
			ar & m_RBC;

			// rigid body displacements
			ar & m_RDC;

			// rigid body initial conditions
			ar & m_RIC;
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidSystem::Activate()
{
	// rigid body displacements
	for (int i=0; i<(int) m_RDC.size(); ++i)
	{
		FERigidPrescribedBC& rc = *m_RDC[i];
		if (rc.IsActive()) rc.Activate();
	}

	// fixed rigid body dofs
	for (int i=0; i<(int) m_RBC.size(); ++i)
	{
		FERigidFixedBC& rc = *m_RBC[i];
		if (rc.IsActive()) rc.Activate();
	}

	// initial rigid conditions
	for (int i=0; i<(int) m_RIC.size(); ++i)
	{
		FERigidIC& ric = *m_RIC[i];
		if (ric.IsActive()) ric.Activate();
	}
}

//-----------------------------------------------------------------------------
bool FERigidSystem::Init()
{
	// create the rigid bodies from the rigid material list
	if (CreateObjects() == false) return false;

	// the rigid body constraints are still associated with the rigid materials
	// so we now associate them with the rigid bodies
	// NOTE: This is now done in the Init() member function for each class
	FEModel& fem = m_fem;
	for (int i=0; i<(int) m_RBC.size(); ++i)
	{
		FERigidFixedBC& BC = *m_RBC[i];
		if (BC.Init() == false) return false;
	}
	for (int i=0; i<(int) m_RDC.size(); ++i)
	{
		FERigidPrescribedBC& DC = *m_RDC[i];
		if (DC.Init() == false) return false;
	}
	for (int i=0; i<(int) m_RIC.size(); ++i)
	{
		FERigidIC& IC = *m_RIC[i];
		if (IC.Init() == false) return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! In FEBio rigid bodies are defined implicitly through a list of rigid materials.
//! However, a rigid body can be composed of multiple rigid materials, so we can't
//! just create a rigid body for each rigid material. We must look at the connectivity
//! of the materials as well.
bool FERigidSystem::CreateObjects()
{
	FEModel& fem = m_fem;
	FEMesh& mesh = fem.GetMesh();

	// count the number of rigid materials
	int NMAT = fem.Materials();
	int nrm = 0;
	for (int i=0; i<NMAT; ++i)
	{
		if (dynamic_cast<FERigidMaterial*>(fem.GetMaterial(i))) nrm++;
	}
	
	// make sure there are rigid materials
	if (nrm == 0) return true;

	// First we need to figure out how many rigid bodies there are.
	// This is not the same as rigid materials, since a rigid body
	// may be composed of different rigid materials (similarly to a deformable
	// body that may contain different materials). Although there can
	// only be one deformable mesh, there can be several rigid bodies.

	// The mrb array will contain an index to the rigid body the material
	// is attached to. For now, we assume one rigid body per rigid material.
	vector<int> mrb(NMAT);
	int n = 0;
	for (int i=0; i<NMAT; ++i)
	{
		if (dynamic_cast<FERigidMaterial*>(fem.GetMaterial(i))) mrb[i] = n++;
		else mrb[i] = -1;
	}

	// Next, we assign to all nodes a rigid node number
	// This number is preliminary since rigid materials can be merged
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		FERigidMaterial* pmat = dynamic_cast<FERigidMaterial*>(dom.GetMaterial());
		if (pmat)
		{
			for (int i=0; i<dom.Elements(); ++i)
			{
				FEElement& el = dom.ElementRef(i);
				for (int j=0; j<el.Nodes(); ++j)
				{
					int n = el.m_node[j];
					FENode& node = mesh.Node(n);
					node.m_rid = pmat->GetID() - 1;
				}
			}
		}
	}

	// now we can merge rigid materials
	// if a rigid element has two nodes that connect to two different
	// rigid materials we need to merge. 
	bool bdone;
	do
	{
		bdone = true;
		for (int nd=0; nd<mesh.Domains(); ++nd)
		{
			FEDomain& dom = mesh.Domain(nd);
			FERigidMaterial* pmat = dynamic_cast<FERigidMaterial*>(dom.GetMaterial());
			if (pmat)
			{
				for (int i=0; i<dom.Elements(); ++i)
				{
					FEElement& el = dom.ElementRef(i);

					int m = mesh.Node(el.m_node[0]).m_rid;
					for (int j=1; j<el.Nodes(); ++j)
					{
						int n = mesh.Node(el.m_node[j]).m_rid;
						if (mrb[n] != mrb[m])
						{
							if (mrb[n]<mrb[m]) mrb[m] = mrb[n]; else mrb[n] = mrb[m];
							bdone = false;
						}
					}
				}
			}
		}
	}
	while (!bdone);

	// since we may have lost a rigid body in the merge process
	// we reindex the RB's.
	vector<int> mrc; mrc.assign(NMAT, -1);
	for (int i=0; i<NMAT; ++i) if (mrb[i] >= 0) mrc[mrb[i]] = 0;
	int nrb = 0;
	for (int i=0; i<NMAT; ++i)
	{
		if (mrc[i] == 0) mrc[i] = nrb++;
	}

	for (int i=0; i<NMAT; ++i) 
	{
		if (mrb[i] >= 0) mrb[i] = mrc[mrb[i]];
	}

	// set rigid body index for materials
	for (int i=0; i<NMAT; ++i)
	{
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(i));
		if (pm)	
		{
			pm->SetRigidBodyID(mrb[i]);
		}
	}

	// assign rigid body index to nodes
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid >= 0) 
		{
			node.m_rid = mrb[ node.m_rid ];
			node.m_ra = node.m_r0;
		}
	}

	// let's clear all rigid bodies
	if (m_RB.empty() == false)
	{
		for (int i=0; i<(int) m_RB.size(); ++i) delete m_RB[i]; 
		m_RB.clear();
	}

	// Ok, we now know how many rigid bodies there are
	// so let's create them
	for (int i=0; i<nrb; ++i)
	{
		// create a new rigid body
		FERigidBody* prb = new FERigidBody(&fem);
		prb->m_nID = i;

		// Since a rigid body may contain several rigid materials
		// we find the first material that this body has and use
		// that materials data to set up the rigid body data
		int j;
		FERigidMaterial* pm = 0;
		for (j=0; j<NMAT; ++j)
		{
			pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(j));
			if (pm && (pm->GetRigidBodyID() == i)) break;
		}
		if (j >= NMAT) return false;
		prb->m_mat = j;
		prb->SetName(pm->GetName());

		// add it to the pile
		m_RB.push_back(prb);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Reset rigid system.
bool FERigidSystem::Reset()
{
	int nrb = (int)m_RB.size();
	for (int i=0; i<nrb; ++i) m_RB[i]->Reset();
	return true;
}

//-----------------------------------------------------------------------------
//! Find a parameter
double* FERigidSystem::FindParameter(int nmat, ParamString& sz, int index)
{
	// the rigid bodies are dealt with differently
	int nrb = (int)m_RB.size();
	for (int i=0; i<nrb; ++i)
	{
		FERigidBody& ob = *m_RB[i];
		if (ob.GetMaterialID() == nmat)
		{
			FEParam* pp = ob.FindParameter(sz);
			if (pp) return pp->pvalue<double>(index);
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
// find a parameter value
FEParamValue FERigidSystem::GetParameterValue(const ParamString& paramString)
{
	ParamString next = paramString.next();
	if (next == "rigidbody")
	{
		int NRB = Objects();

		if (next.Index() >= 0)
		{
			int index = next.Index();
			if ((index >= 0) && (index < NRB))
			{
				ParamString paramName = next.next();

				FERigidBody* ob = Object(index);
				FEParam* pp = ob->FindParameter(paramName);
				if (pp) return GetParameterComponent(paramName.last(), pp);
			}
		}
		else
		{
			FEMaterial* mat = 0;
			if (next.IDString()) mat = m_fem.FindMaterial(next.IDString());
			if ((mat != 0) && (dynamic_cast<FERigidMaterial*>(mat)))
			{
				ParamString paramName = next.next();

				// the rigid bodies are dealt with differently
				int nmat = mat->GetID() - 1;
				for (int i = 0; i < NRB; ++i)
				{
					FERigidBody* ob = Object(i);
					if (ob && (ob->GetMaterialID() == nmat))
					{
						FEParam* pp = ob->FindParameter(paramName);
						if (pp) return GetParameterComponent(paramName.last(), pp);
					}
				}
			}
		}
	}

	return FEParamValue();
}

//-----------------------------------------------------------------------------
//! Call this function after the rigid body kinematics are updated to ensure
//! that the mesh (or at least the part of the mesh corresponding to the rigid
//! bodies) is updated.
void FERigidSystem::UpdateMesh()
{
	FEMesh& mesh = m_fem.GetMesh();
	int NRB = Objects();
	for (int i=0; i<NRB; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *Object(i);

		// update the mesh' nodes
		int N = mesh.Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = mesh.Node(i);
			if (node.m_rid == RB.m_nID)
			{
				vec3d a0 = node.m_ra - RB.m_r0;
				vec3d at = RB.GetRotation()*a0;
				node.m_rt = RB.m_rt + at;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidSystem::BuildMatrixProfile(FEGlobalMatrix& G)
{
	if (Objects())
	{
		vector<int> lm(6);
		int nrb = Objects();
		for (int i=0; i<nrb; ++i)
		{
			FERigidBody& rb = *Object(i);
			for (int j=0; j<6; ++j) lm[j] = rb.m_LM[j];
			G.build_add(lm);
		}
	}
}
