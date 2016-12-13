#include "stdafx.h"
#include "FEPeriodicLinearConstraint2O.h"
#include <FECore/FEModel.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FELinearConstraint.h>
#include <FECore/FEMesh.h>
#include <FECore/FESurface.h>

FEPeriodicLinearConstraint2O::NodeSetSet::NodeSetSet()
{
}

FEPeriodicLinearConstraint2O::NodeSetSet::NodeSetSet(const FEPeriodicLinearConstraint2O::NodeSetSet& nss)
{
	master = nss.master;
	slave = nss.slave;
}

void FEPeriodicLinearConstraint2O::NodeSetSet::operator = (const FEPeriodicLinearConstraint2O::NodeSetSet& nss)
{
	master = nss.master;
	slave = nss.slave;
}


FEPeriodicLinearConstraint2O::FEPeriodicLinearConstraint2O()
{
}

FEPeriodicLinearConstraint2O::~FEPeriodicLinearConstraint2O()
{
}

void FEPeriodicLinearConstraint2O::AddNodeSetPair(const FENodeSet& ms, const FENodeSet& ss, bool push_back)
{
	NodeSetSet sp;
	sp.master = ms;
	sp.slave.push_back(ss);
	if (push_back) m_set.push_back(sp); else m_set.insert(m_set.begin(), sp);
}

void FEPeriodicLinearConstraint2O::AddNodeSetSet(FENodeSet** set, int count, bool push_back)
{
	if (count <= 1) return;
	NodeSetSet sp;
	sp.master = *set[0];
	for (int i=1; i<count; ++i)sp.slave.push_back(*set[i]);

	if (push_back) m_set.push_back(sp); else m_set.insert(m_set.begin(), sp);
}

bool FEPeriodicLinearConstraint2O::GenerateConstraints(FEModel* fem)
{
	// get the model's mesh
	FEMesh& mesh = fem->GetMesh();

	// make sure there is a list of sets
	if (m_set.empty()) return true;

	// tag all the mesh' nodes
	// 0 = no constraint applied yet
	// 1 = already constrained (or cannot be constrained)
	vector<int> tag(mesh.Nodes(), 0);

	// tag the exlude list
	FENodeSet& ns = m_exclude;
	for (int i = 0; i<ns.size(); ++i) tag[ns[i]] = 1;

	// get the linear constraint manager
	FELinearConstraintManager& LCM = fem->GetLinearConstraintManager();

	// loop over all the sets
	for (int n = 0; n<(int)m_set.size(); ++n)
	{
		// extract the nodes from the surfaces
		FENodeSet& master = m_set[n].master;
		vector<FENodeSet>& slaveSet = m_set[n].slave;

		// loop over all slave surfaces
		for (int m=0; m<(int)slaveSet.size(); ++m)
		{
			FENodeSet& slave = slaveSet[m];

			// first, we find for each master node, its closest slave node
			vector<pair<int, int> > nodes;

			// loop over all slave nodes
			for (int i = 0; i<slave.size(); ++i)
			{
				// get the slave node position
				vec3d& rm = slave.Node(i)->m_r0;

				// find the closest master node
				int m = -1;
				double Dmin = 0.0;
				for (int j = 0; j<master.size(); ++j)
				{
					vec3d& rs = master.Node(j)->m_r0;
					double D = (rm - rs)*(rm - rs);
					if ((D < Dmin) || (m == -1))
					{
						Dmin = D;
						m = j;
					}
				}
				if (m == -1) return false;

				pair<int, int> p(slave[i], master[m]);
				nodes.push_back(p);
			}

			// setup all the linear constraints
			for (int i = 0; i<(int)nodes.size(); ++i)
			{
				pair<int, int> p = nodes[i];

				if (p.first == p.second)
				{
					return false;
				}

				// make sure the slave is not excluded 
				if (tag[p.first] == 0)
				{
					// do one constraint for x, y, z
					for (int j = 0; j<3; ++j)
					{
						FELinearConstraint lc(fem);
						lc.SetMasterDOF(j, p.first);

						lc.AddSlaveDof(j, p.second, 1.0);

						LCM.AddLinearConstraint(lc);
					}

					// don't forget to tag the node
					tag[p.first] = 1;
				}
			}
		}
	}

	return true;
}
