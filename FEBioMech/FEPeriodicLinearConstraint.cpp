#include "stdafx.h"
#include "FEPeriodicLinearConstraint.h"
#include <FECore/FEModel.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FELinearConstraint.h>
#include <FECore/FEMesh.h>
#include <FECore/FESurface.h>

FEPeriodicLinearConstraint::FEPeriodicLinearConstraint()
{
	m_refNode = -1;
}

FEPeriodicLinearConstraint::~FEPeriodicLinearConstraint()
{
}

void FEPeriodicLinearConstraint::AddNodeSetPair(const FENodeSet& ms, const FENodeSet& ss, bool push_back)
{
	NodeSetPair sp;
	sp.master = ms;
	sp.slave = ss;
	if (push_back) m_set.push_back(sp); else m_set.insert(m_set.begin(), sp);
}

bool FEPeriodicLinearConstraint::GenerateConstraints(FEModel* fem)
{
	// get the model's mesh
	FEMesh& mesh = fem->GetMesh();

	// make sure there is a list of sets
	if (m_set.empty()) return true;

	// check the reference node
	if ((m_refNode < 0) || (m_refNode >= mesh.Nodes())) return false;

	// tag all the mesh' nodes
	// 0 = no constraint applied yet
	// 1 = already constrained (or cannot be constrained)
	vector<int> tag(mesh.Nodes(), 0);

	// tag the exlude list
	FENodeSet& ns = m_exclude;
	for (int i=0; i<ns.size(); ++i) tag[ ns[i] ] = 1;

	// get the linear constraint manager
	FELinearConstraintManager& LCM = fem->GetLinearConstraintManager();

	// loop over all the pairs
	for (int n=0; n<(int)m_set.size(); ++n)
	{
		// extract the nodes from the surfaces
		FENodeSet& master = m_set[n].master;
		FENodeSet& slave  = m_set[n].slave;

		// first, we find for each master node, it's closest slave node
		vector<pair<int, int> > nodes;

		// loop over all master nodes
		// (These are the dofs that will be eliminated)
		for (int i=0; i<master.size(); ++i)
		{
			// get the master node position
			vec3d& rm = master.Node(i)->m_r0;

			// find the closest slave node
			int m = -1;
			double Dmin = 0.0;
			for (int j=0; j<slave.size(); ++j)
			{
				vec3d& rs = slave.Node(j)->m_r0;
				double D = (rm - rs)*(rm - rs);
				if ((D < Dmin) || (m == -1))
				{
					Dmin = D;
					m = j;
				}
			}

			pair<int, int> p(master[i], slave[m]);
			nodes.push_back(p);
		}

		// find the index of the reference node in the slave node list
		int mRef = -1;
		for (int i=0; i<(int) nodes.size(); ++i)
		{
			if (nodes[i].second == m_refNode)
			{
				mRef = nodes[i].first;
				break;
			}
		}
		// make sure the reference node was found in the slave list
		if (mRef == -1) return false;

		// setup all the linear constraints
		for (int i=0; i<(int) nodes.size(); ++i)
		{
			pair<int,int> p = nodes[i];

			// make sure the master is not excluded 
			if (tag[p.first] == 0)
			{
				// do one constraint for x, y, z
				for (int j=0; j<3; ++j)
				{
					FELinearConstraint lc(fem);
					lc.SetMasterDOF(j, p.first);

					lc.AddSlaveDof(j, p.second ,  1.0);
					lc.AddSlaveDof(j, mRef     ,  1.0);
					lc.AddSlaveDof(j, m_refNode, -1.0);

					LCM.AddLinearConstraint(lc);
				}

				// don't forget to tag the node
				tag[p.first] = 1;
			}
		}
	}

	return true;
}
