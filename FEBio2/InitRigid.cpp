#include "stdafx.h"
#include "fem.h"
#include "FEBioLib/FERigid.h"
#include "FEBioLib/FERigidJoint.h"
#include "FEBioLib/FEElasticShellDomain.h"
#include "FEBioLib/FEElasticSolidDomain.h"
#include "FEBioLib/log.h"

//-----------------------------------------------------------------------------
//! This function creates the rigid bodies by analyzing the rigid materials
//! and the mesh in the model. 
//!
bool FEM::CreateRigidBodies()
{
	int i, j, n, m, nd;
	// count the number of rigid materials
	int nrm = 0;
	for (i=0; i<Materials(); ++i)
	{
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(i));
		if (pm) nrm++;
	}
	
	// make sure there are rigid materials
	if (nrm == 0) return true;

	// First we need to figure out how many rigid bodies there are.
	// This is not the same as rigid materials, since a rigid body
	// may be composed of different rigid materials (similarly to a deformable
	// body that may contain different materials). Although there can
	// only be one deformable mesh, there can be several rigid bodies.

	// The mrb array will contain an index to the rigid body the material
	// is attached too.
	vector<int> mrb(Materials());
	n = 0;
	for (i=0; i<Materials(); ++i)
	{
		if (dynamic_cast<FERigidMaterial*>(GetMaterial(i))) mrb[i] = n++;
		else mrb[i] = -1;
	}

	// Next, we assign to all nodes a rigid node number
	// This number is preliminary since rigid materials can be merged
	for (nd = 0; nd < m_mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(el.GetMatID()));
				if (pm)
				{
					el.m_nrigid = el.GetMatID();
					for (j=0; j<el.Nodes(); ++j)
					{
						n = el.m_node[j];
						FENode& node = m_mesh.Node(n);
						node.m_rid = el.GetMatID();
					}
				}
				else el.m_nrigid = -1;
			}
		}

		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(el.GetMatID()));
				if (pm)
				{
					el.m_nrigid = el.GetMatID();
					for (j=0; j<el.Nodes(); ++j)
					{
						n = el.m_node[j];
						FENode& node = m_mesh.Node(n);
						node.m_rid = el.GetMatID();
					}		
				}
				else el.m_nrigid = -1;
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
		for (nd=0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if (el.m_nrigid >= 0)
					{
						m = m_mesh.Node(el.m_node[0]).m_rid;
						for (j=1; j<el.Nodes(); ++j)
						{
							n = m_mesh.Node(el.m_node[j]).m_rid;
							if (mrb[n] != mrb[m])
							{
								if (mrb[n]<mrb[m]) mrb[m] = mrb[n]; else mrb[n] = mrb[m];
								bdone = false;
							}
						}
					}
				}
			}

			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if (el.m_nrigid >= 0)
					{
						m = m_mesh.Node(el.m_node[0]).m_rid;
						for (j=1; j<el.Nodes(); ++j)
						{
							n = m_mesh.Node(el.m_node[j]).m_rid;
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
	}
	while (!bdone);

	// since we may have lost a rigid body in the merge process
	// we reindex the RB's.
	int nmat = Materials();
	vector<int> mrc; mrc.assign(nmat, -1);
	for (i=0; i<nmat; ++i) if (mrb[i] >= 0) mrc[mrb[i]] = 0;
	int nrb = 0;
	for (i=0; i<nmat; ++i)
	{
		if (mrc[i] == 0) mrc[i] = nrb++;
	}

	for (i=0; i<nmat; ++i) 
	{
		if (mrb[i] >= 0) mrb[i] = mrc[mrb[i]];
	}

	// set rigid body index for materials
	for (i=0; i<Materials(); ++i)
	{
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*> (GetMaterial(i));
		if (pm)	
		{
			pm->m_nRB = mrb[i];
		}
	}

	// assign rigid body index to rigid elements
	for (nd=0; nd<m_mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*> (GetMaterial(el.GetMatID()));
				if (pm)
					el.m_nrigid = pm->m_nRB;
				else
					el.m_nrigid = -1;
			}
		}

		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				FERigidMaterial* pm = dynamic_cast<FERigidMaterial*> (GetMaterial(el.GetMatID()));
				if (pm)
					el.m_nrigid = pm->m_nRB;
				else
					el.m_nrigid = -1;
			}
		}
	}

	// assign rigid body index to nodes
	for (i=0; i<m_mesh.Nodes(); ++i)
	{
		FENode& node = m_mesh.Node(i);
		if (node.m_rid >= 0) node.m_rid = mrb[ node.m_rid ];
	}

	// Ok, we now know how many rigid bodies there are
	// so let's create them
	m_Obj.clear();
	for (i=0; i<nrb; ++i)
	{
		// create a new rigid body
		FERigidBody* prb = new FERigidBody(this);
		prb->m_nID = i;

		// Since a rigid body may contain several rigid materials
		// we find the first material that this body has and use
		// that materials data to set up the rigid body data
		FERigidMaterial* pm = 0;
		for (j=0; j<Materials(); ++j)
		{
			pm = dynamic_cast<FERigidMaterial*> (GetMaterial(j));

			if (pm && (pm->m_nRB == i))	break;
		}
		assert(j<Materials());
		prb->m_mat = j;

		// initialize center of mass
		if (pm->m_com == 1)
		{
			// grab the com from the material
			prb->m_r0 = prb->m_rt = pm->m_rc;
		}
		else
		{
			// calculate the com
			prb->UpdateCOM();
		}

		// add it to the pile
		m_Obj.push_back(prb);
	}

	// set up rigid joints
	if (!m_NLC.empty())
	{
		FERigidMaterial* pm;
		int NC = m_NLC.size();
		for (i=0; i<NC; ++i)
		{
			FENLConstraint* plc = m_NLC[i];
			if (plc->Type() == FE_RIGID_JOINT)
			{
				FERigidJoint& rj = dynamic_cast<FERigidJoint&>(*plc);
				rj.m_F = vec3d(0,0,0);

				pm = dynamic_cast<FERigidMaterial*> (GetMaterial(rj.m_nRBa));
				if (pm == 0)
				{
					clog.printbox("FATAL ERROR", "Rigid joint %d does not connect two rigid bodies\n", i+1);
					return false;
				}
				rj.m_nRBa = pm->m_nRB;

				pm = dynamic_cast<FERigidMaterial*> (GetMaterial(rj.m_nRBb));
				if (pm == 0)
				{
					clog.printbox("FATAL ERROR", "Rigid joint %d does not connect two rigid bodies\n", i+1);
					return false;
				}
				rj.m_nRBb = pm->m_nRB;

				FERigidBody& ra = dynamic_cast<FERigidBody&>(*m_Obj[rj.m_nRBa]);
				FERigidBody& rb = dynamic_cast<FERigidBody&>(*m_Obj[rj.m_nRBb]);

				rj.m_qa0 = rj.m_q0 - ra.m_r0;
				rj.m_qb0 = rj.m_q0 - rb.m_r0;
			}
		}
	}

	// overwrite rigid nodes degrees of freedom
	// We do this so that these dofs do not
	// get equation numbers assigned to them. Later we'll assign
	// the rigid dofs equations numbers to these nodes
	for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_bshell = false;
	for (nd = 0; nd<m_mesh.Domains(); ++nd)
	{
		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				if (el.m_nrigid < 0)
				{
					int n = el.Nodes();
					for (j=0; j<n; ++j) m_mesh.Node(el.m_node[j]).m_bshell = true;
				}
			}
		}
	}

	// The following fixes the degrees of freedom for rigid nodes.
	// Note that also the rotational degrees of freedom are fixed
	// for rigid nodes that do not belong to a non-rigid shell element.
	for (i=0; i<m_mesh.Nodes(); ++i)
	{
		FENode& node = m_mesh.Node(i);
		if (node.m_rid >= 0)
		{
			node.m_ID[DOF_X] = -1;
			node.m_ID[DOF_Y] = -1;
			node.m_ID[DOF_Z] = -1;
			if (node.m_bshell == false)
			{
				node.m_ID[DOF_U] = -1;
				node.m_ID[DOF_V] = -1;
				node.m_ID[DOF_W] = -1;
			}
		}
	}

	// assign correct rigid body ID's to rigid nodes
	for (i=0; i<(int) m_RN.size(); ++i)
	{
		FERigidNode& rn = *m_RN[i];
		rn.rid = mrb[rn.rid];
	}

	// let's find all rigid surface elements
	// a surface element is rigid when it has no free nodes
	for (int is = 0; is < (int) m_SL.size(); ++is)
	{
		FESurfaceLoad* ps = m_SL[is];
		for (i=0; i<ps->Surface().Elements(); ++i)
		{
			FESurfaceElement& el = ps->Surface().Element(i);
			int N = el.Nodes();
			el.m_nrigid = 0;
			for (j=0; j<N; ++j) 
			{
				FENode& node = m_mesh.Node(el.m_node[j]);
				if (node.m_rid < 0) 
				{
					el.m_nrigid = -1;
					break;
				}
			}
		}
	}

	// the rigid body constraints are still associated with the rigid materials
	// so we now associate them with the rigid bodies
	for (i=0; i<(int) m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = *m_RDC[i];
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(DC.id-1));
		DC.id = pm->m_nRB;
	}
	for (i=0; i<(int) m_RFC.size(); ++i)
	{
		FERigidBodyForce& FC = *m_RFC[i];
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(GetMaterial(FC.id-1));
		FC.id = pm->m_nRB;
	}

	// set the rigid body parents
	for (i=0; i<(int) m_Obj.size(); ++i)
	{
		FERigidBody& rb = dynamic_cast<FERigidBody&>(*m_Obj[i]);
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_MAT[rb.m_mat]);
		assert(pm);
		if (pm->m_pmid > -1)
		{
			FERigidMaterial* ppm = dynamic_cast<FERigidMaterial*>(m_MAT[pm->m_pmid-1]);
			assert(ppm);
			FERigidBody& prb = dynamic_cast<FERigidBody&>(*m_Obj[ppm->m_nRB]);
			rb.m_prb = &prb;

			// we also need to open up all the RB's degree of freedoms
			pm->m_bc[0] = 1;
			pm->m_bc[1] = 1;
			pm->m_bc[2] = 1;
			pm->m_bc[3] = 1;
			pm->m_bc[4] = 1;
			pm->m_bc[5] = 1;
		}
	}

	return true;
}
