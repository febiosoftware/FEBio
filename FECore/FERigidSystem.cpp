#include "stdafx.h"
#include "FERigidSystem.h"
#include "FEModel.h"

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
int FERigidSystem::Objects() const
{
	return (int) m_RB.size();
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
	int i;
	for (i=0; i<m_RDC.size(); ++i) delete m_RDC[i]; m_RDC.clear();
	for (i=0; i<m_RBV.size(); ++i) delete m_RBV[i]; m_RBV.clear();
	for (i=0; i<m_RBW.size(); ++i) delete m_RBW[i]; m_RBW.clear();
	for (i=0; i<m_RN.size (); ++i) delete m_RN [i]; m_RN.clear ();
}

//-----------------------------------------------------------------------------
void FERigidSystem::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		// rigid nodes
		ar << (int) m_RN.size();
		for (int i=0; i<(int) m_RN.size(); ++i)
		{
			FERigidNode& rn = *m_RN[i];
			rn.Serialize(ar);
		}

		// fixed rigid body dofs
		ar << (int) m_RBC.size();
		for (int i=0; i<(int) m_RBC.size(); ++i)
		{
			FERigidBodyFixedBC& bc = *m_RBC[i];
			bc.Serialize(ar);
		}

		// rigid body displacements
		ar << (int) m_RDC.size();
		for (int i=0; i<(int) m_RDC.size(); ++i)
		{
			FERigidBodyDisplacement& dc = *m_RDC[i];
			dc.Serialize(ar);
		}
	}
	else
	{
		// rigid nodes
		int n = 0;
		ar >> n;
		m_RN.clear();
		for (int i=0; i<n; ++i)
		{
			FERigidNode* prn = new FERigidNode(&m_fem);
			prn->Serialize(ar);
			if (prn->IsActive()) prn->Activate(); else prn->Deactivate();
			m_RN.push_back(prn);
		}

		// fixed rigid body dofs
		ar >> n;
		m_RBC.clear();
		for (int i=0; i<n; ++i)
		{
			FERigidBodyFixedBC* pbc = new FERigidBodyFixedBC(&m_fem);
			pbc->Serialize(ar);
			if (pbc->IsActive()) pbc->Activate(); else pbc->Deactivate();
			m_RBC.push_back(pbc);
		}

		// rigid body displacements
		ar >> n;
		m_RDC.clear();
		for (int i=0; i<n; ++i)
		{
			FERigidBodyDisplacement* pdc = new FERigidBodyDisplacement(&m_fem);
			pdc->Serialize(ar);
			if (pdc->IsActive()) pdc->Activate(); else pdc->Deactivate();
			m_RDC.push_back(pdc);
		}
	}
}

//-----------------------------------------------------------------------------
//! Find a BC based on its ID. This is needed for restarts.
FEModelComponent* FERigidSystem::FindModelComponent(int nid)
{
	int i;
	for (i=0; i<(int) m_RBC.size(); ++i) if (m_RBC[i]->GetClassID() == nid) return m_RBC[i];
	for (i=0; i<(int) m_RDC.size(); ++i) if (m_RDC[i]->GetClassID() == nid) return m_RDC[i];
	for (i=0; i<(int) m_RN.size (); ++i) if (m_RN [i]->GetClassID() == nid) return m_RN [i];
	return 0;
}

//-----------------------------------------------------------------------------
void FERigidSystem::Activate()
{
	// rigid nodes
	for (int i=0; i<(int) m_RN.size(); ++i)
	{
		FERigidNode& rn = *m_RN[i];
		if (rn.IsActive()) rn.Activate();
	}

	// rigid body displacements
	for (int i=0; i<(int) m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& rc = *m_RDC[i];
		if (rc.IsActive()) rc.Activate();
	}

	// fixed rigid body dofs
	for (int i=0; i<(int) m_RBC.size(); ++i)
	{
		FERigidBodyFixedBC& rc = *m_RBC[i];
		if (rc.IsActive()) rc.Activate();
	}

	// initial rigid velocity
	for (int i=0; i<(int) m_RBV.size(); ++i)
	{
		FERigidBodyVelocity& RV = *m_RBV[i];
		if (RV.IsActive()) RV.Activate();
	}

	// initial rigid angular velocity
	for (int i=0; i<(int) m_RBW.size(); ++i)
	{
		FERigidBodyAngularVelocity& RW = *m_RBW[i];
		if (RW.IsActive()) RW.Activate();
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
		FERigidBodyFixedBC& BC = *m_RBC[i];
		if (BC.Init() == false) return false;
	}
	for (int i=0; i<(int) m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *m_RDC[i];
		if (DC.Init() == false) return false;
	}
	for (int i=0; i<(int) m_RBV.size(); ++i)
	{
		FERigidBodyVelocity& RV = *m_RBV[i];
		if (RV.Init() == false) return false;
	}
	for (int i=0; i<(int) m_RBW.size(); ++i)
	{
		FERigidBodyAngularVelocity& RW = *m_RBW[i];
		if (RW.Init() == false) return false;
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

	// Find the nodes that are on a non-rigid shell. 
	// These nodes will be assigned rotational degrees of freedom
	// TODO: I'm not sure if this is a good place to do this.
	for (int i=0; i<mesh.Nodes(); ++i) mesh.Node(i).m_bshell = false;
	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		if (dom.Class() == FE_DOMAIN_SHELL)
		{
			FEMaterial* pmat = dom.GetMaterial();
			if (pmat->IsRigid() == false)
			{
				int N = dom.Elements();
				for (int i=0; i<N; ++i)
				{
					FEElement& el = dom.ElementRef(i);
					int n = el.Nodes();
					for (int j=0; j<n; ++j) mesh.Node(el.m_node[j]).m_bshell = true;
				}
			}
		}
	}

	// count the number of rigid materials
	int NMAT = fem.Materials();
	int nrm = 0;
	for (int i=0; i<NMAT; ++i)
	{
		if (fem.GetMaterial(i)->IsRigid()) nrm++;
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
		if (fem.GetMaterial(i)->IsRigid()) mrb[i] = n++;
		else mrb[i] = -1;
	}

	// Next, we assign to all nodes a rigid node number
	// This number is preliminary since rigid materials can be merged
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		FEMaterial* pmat = dom.GetMaterial();
		if (pmat->IsRigid())
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
			FEMaterial* pmat = dom.GetMaterial();
			if (pmat->IsRigid())
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
		FEMaterial* pm = fem.GetMaterial(i);
		if (pm->IsRigid())	
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
		FEMaterial* pm = 0;
		for (j=0; j<NMAT; ++j)
		{
			pm = fem.GetMaterial(j);
			if (pm && (pm->GetRigidBodyID() == i))	break;
		}
		if (j >= NMAT) return false;
		prb->m_mat = j;

		// add it to the pile
		m_RB.push_back(prb);
	}

	// assign correct rigid body ID's to rigid nodes
	for (int i=0; i<(int) m_RN.size(); ++i)
	{
		FERigidNode& rn = *m_RN[i];
		rn.rid = mrb[rn.rid];
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Reset rigid system.
bool FERigidSystem::Reset()
{
	int nrb = m_RB.size();
	for (int i=0; i<nrb; ++i) m_RB[i]->Reset();
	return true;
}

//-----------------------------------------------------------------------------
//! Copy some data to the dump stream for running restarts
void FERigidSystem::ShallowCopy(DumpStream& dmp, bool bsave)
{
	for (int i=0; i<(int) m_RB.size(); ++i) m_RB[i]->ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
//! Find a parameter
double* FERigidSystem::FindParameter(int nmat, ParamString& sz, int index)
{
	// the rigid bodies are dealt with differently
	int nrb = m_RB.size();
	for (int i=0; i<nrb; ++i)
	{
		FEObject& ob = *m_RB[i];
		if (ob.GetMaterialID() == nmat)
		{
			FEParam* pp = ob.GetParameter(sz);
			if (pp) return pp->pvalue<double>(index);
		}
	}

	return 0;
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
				vec3d a0 = node.m_r0 - RB.m_r0;
				vec3d at = RB.m_qt*a0;
				node.m_rt = RB.m_rt + at;
			}
		}
	}
}
