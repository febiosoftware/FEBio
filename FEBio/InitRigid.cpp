#include "stdafx.h"
#include "fem.h"
#include "FERigid.h"
#include "stack.h"
#include "log.h"

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEM::InitRigidBodies
//  Initializes rigid body data
//

bool FEM::InitRigidBodies()
{
	int i, j, n, m, nd;
	// count the number of rigid materials
	m_nrm = 0;
	for (i=0; i<Materials(); ++i)
	{
		FERigid* pm = dynamic_cast<FERigid*>(GetMaterial(i));
		if (pm) m_nrm++;
	}
	
	// make sure there are rigid materials
	if (m_nrm == 0) return true;

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
		if (dynamic_cast<FERigid*>(GetMaterial(i))) mrb[i] = n++;
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
				FERigid* pm = dynamic_cast<FERigid*>(GetMaterial(el.GetMatID()));
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
				FERigid* pm = dynamic_cast<FERigid*>(GetMaterial(el.GetMatID()));
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

	// now let's count how many rigid bodies we have
	m_nrb = -1;
	for (i=0; i<Materials(); ++i) if (mrb[i] > m_nrb) m_nrb = mrb[i];
	m_nrb++;

	// set rigid body index for materials
	for (i=0; i<Materials(); ++i)
	{
		FERigid* pm = dynamic_cast<FERigid*> (GetMaterial(i));
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
				FERigid* pm = dynamic_cast<FERigid*> (GetMaterial(el.GetMatID()));
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
				FERigid* pm = dynamic_cast<FERigid*> (GetMaterial(el.GetMatID()));
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
	m_RB.resize(m_nrb);
	for (i=0; i<m_nrb; ++i)
	{
		// attach the rigid body do this FEM
		m_RB[i].AttachToFEM(this);

		// zero total displacements
		m_RB[i].m_Ut[0] = m_RB[i].m_Up[0] = 0;
		m_RB[i].m_Ut[1] = m_RB[i].m_Up[1] = 0;
		m_RB[i].m_Ut[2] = m_RB[i].m_Up[2] = 0;
		m_RB[i].m_Ut[3] = m_RB[i].m_Up[3] = 0;
		m_RB[i].m_Ut[4] = m_RB[i].m_Up[4] = 0;
		m_RB[i].m_Ut[5] = m_RB[i].m_Up[5] = 0;

		m_RB[i].m_nID = i;

		// initialize orientation
		m_RB[i].m_qt = quatd(0, vec3d(0,0,1));

		// Since a rigid body may contain several rigid materials
		// we find the first material that this body has and use
		// that materials data to set up the rigid body data
		FERigid* pm = 0;
		for (j=0; j<Materials(); ++j)
		{
			pm = dynamic_cast<FERigid*> (GetMaterial(j));

			if (pm && (pm->m_nRB == i))	break;
		}

		m_RB[i].m_mat = j;

/*		// initialize constraints
		for (j=0; j<6; ++j)
		{
			m_RB[i].m_bc[j] = pm->m_bc[j];
		}
*/
		// initialize center of mass
		if (pm->m_com == 1)
		{
			// grab the com from the material
			m_RB[i].m_r0 = m_RB[i].m_rt = pm->m_rc;
		}
		else
		{
			// calculate the com
			m_RB[i].Update();
		}
	}

	// get the logfile
	Logfile& log = GetLogfile();

	// set up rigid joints
	if (m_nrj > 0)
	{
		FERigid* pm;
		for (i=0; i<m_nrj; ++i)
		{
			FERigidJoint& rj = m_RJ[i];
			rj.m_F = vec3d(0,0,0);

			pm = dynamic_cast<FERigid*> (GetMaterial(rj.m_nRBa));
			if (pm == 0)
			{
				log.printbox("FATAL ERROR", "Rigid joint %d does not connect two rigid bodies\n", i+1);
				return false;
			}
			rj.m_nRBa = pm->m_nRB;

			pm = dynamic_cast<FERigid*> (GetMaterial(rj.m_nRBb));
			if (pm == 0)
			{
				log.printbox("FATAL ERROR", "Rigid joint %d does not connect two rigid bodies\n", i+1);
				return false;
			}
			rj.m_nRBb = pm->m_nRB;

			rj.m_qa0 = rj.m_q0 - m_RB[ rj.m_nRBa ].m_r0;
			rj.m_qb0 = rj.m_q0 - m_RB[ rj.m_nRBb ].m_r0;
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
			node.m_ID[0] = -1;
			node.m_ID[1] = -1;
			node.m_ID[2] = -1;
			if (node.m_bshell == false)
			{
				node.m_ID[3] = -1;
				node.m_ID[4] = -1;
				node.m_ID[5] = -1;
			}
		}
	}

	// assign correct rigid body ID's to rigid nodes
	for (i=0; i<m_RN.size(); ++i)
	{
		FERigidNode& rn = m_RN[i];
		rn.rid = mrb[rn.rid];
	}

	// let's find all rigid surface elements
	// a surface element is rigid when it has no free nodes
	if (m_psurf)
	{
		for (i=0; i<m_psurf->Elements(); ++i)
		{
			FESurfaceElement& el = m_psurf->Element(i);
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
	for (i=0; i<m_RDC.size(); ++i)
	{
		FERigidBodyDisplacement& DC = m_RDC[i];
		FERigid* pm = dynamic_cast<FERigid*>(GetMaterial(DC.id-1));
		DC.id = pm->m_nRB;
	}
	for (i=0; i<m_RFC.size(); ++i)
	{
		FERigidBodyForce& FC = m_RFC[i];
		FERigid* pm = dynamic_cast<FERigid*>(GetMaterial(FC.id-1));
		FC.id = pm->m_nRB;
	}

	return true;
}
