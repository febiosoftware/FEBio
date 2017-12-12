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

// helper function for finding the closest node
int closestNode(FEMesh& mesh, const FENodeSet& set, const vec3d& r);

// helper function for adding the linear constraints
void addLinearConstraint(FEModel& fem, int master, int slave, int nodeA, int nodeB);

// This function generates a linear constraint set based on the definition
// of three surface pairs.
bool FEPeriodicLinearConstraint::GenerateConstraints(FEModel* fem)
{
	// get the model's mesh
	FEMesh& mesh = fem->GetMesh();

	// make sure there is a list of sets
	if (m_set.empty()) return true;

	// make sure we have three sets
	if (m_set.size() != 3) return false;

	// tag all nodes to identify what they are (edge, face, corner)
	int N = mesh.Nodes();
	vector<int> tag(N, 0);

	for (size_t i=0; i<m_set.size(); ++i)
	{
		FENodeSet& ms = m_set[i].master;
		FENodeSet& ss = m_set[i].slave;

		for (int j=0; j<(int)ms.size(); ++j) tag[ms[j]]++;
		for (int j=0; j<(int)ss.size(); ++j) tag[ss[j]]++;
	}

	// flip signs on slave
	for (size_t i = 0; i<m_set.size(); ++i)
	{
		FENodeSet& ss = m_set[i].slave;
		for (int j = 0; j<(int)ss.size(); ++j) 
		{
			int ntag = tag[ss[j]];
			if (ntag > 0) tag[ss[j]] = -ntag;
		}
	}

	// At this point, the following should hold
	// slave nodes: tag < 0, master nodes: tag > 0, interior nodes: tag = 0
	// only one master node should have a value of 3. We make this the reference node
	int refNode = -1;
	for (int i=0; i<N; ++i)
	{
		if (tag[i] == 3)
		{
			assert(refNode == -1);
			refNode = i;
		}
	}
	assert(refNode != -1);
	if (refNode == -1) return false;

	// get the position of the reference node
	vec3d rm = mesh.Node(refNode).m_r0;

	// create the linear constraints for the surface nodes that don't belong to an edge (i.e. tag = 1)
	for (size_t n = 0; n<m_set.size(); ++n)
	{
		FENodeSet& ms = m_set[n].master;
		FENodeSet& ss = m_set[n].slave;

		// find the corresponding reference node on the slave surface
		int mref = closestNode(mesh, ss, rm);

		// make sure this is a corner node
		assert(tag[ss[mref]] == -3);

		// repeat for all slave nodes
		for (int i=0; i<(int)ss.size(); ++i)
		{
			assert(tag[ss[i]] < 0);
			if (tag[ss[i]] == -1)
			{
				// get the slave node position
				vec3d& rs = ss.Node(i)->m_r0;

				// find the closest master node
				int m = closestNode(mesh, ms, rs);
				assert(tag[ms[m]] == 1);

				// add the linear constraints
				addLinearConstraint(*fem, ss[i], ms[m], ss[mref], refNode);
			}
		}
	}

	// extract all 12 edges
	vector<FENodeSet*> surf;
	for (int i=0; i<(int)m_set.size(); ++i)
	{
		surf.push_back(&m_set[i].master);
		surf.push_back(&m_set[i].slave);
	}

	vector<FENodeSet> masterEdges;
	vector<FENodeSet> slaveEdges;
	for (int i=0; i<surf.size(); ++i)
	{
		FENodeSet& s0 = *surf[i];
		for (int j=i+1; j<surf.size(); ++j)
		{
			FENodeSet& s1 = *surf[j];
			vector<int> tmp(N, 0);
			for (int k=0; k<s0.size(); ++k) tmp[s0[k]]++;
			for (int k=0; k<s1.size(); ++k) tmp[s1[k]]++;

			FENodeSet edge(&mesh);
			for (int k=0; k<N; ++k)
			{
				if (tmp[k] == 2) edge.add(k);
			}

			if (edge.size() != 0)
			{
				// see if this is a master edge or not
				// we assume it's a master edge if it connects to the refnode
				bool bmaster = false;
				for (int k=0; k<edge.size(); ++k)
				{
					if (edge[k] == refNode)
					{
						bmaster = true;
						break;
					}
				}

				if (bmaster)
					masterEdges.push_back(edge);
				else
					slaveEdges.push_back(edge);
			}
		}
	}

	// since it is assumed the geometry is a cube, the following must hold
	assert(masterEdges.size() == 3);
	assert(slaveEdges.size() == 9);

	// find the master edge vectors
	vec3d Em[3];
	for (int i=0; i<3; ++i)
	{
		FENodeSet& edge = masterEdges[i];

		// get the edge vector
		Em[i] = edge.Node(0)->m_r0 - edge.Node(1)->m_r0; assert(edge[0] != edge[1]); 
		Em[i].unit();
	}

	// setup linear constraints for edges
	for (int n = 0; n<slaveEdges.size(); ++n)
	{
		FENodeSet& edge = slaveEdges[n];

		// get the edge vector
		vec3d E = edge.Node(0)->m_r0 - edge.Node(1)->m_r0; assert(edge[0] != edge[1]); E.unit();

		// find the corresponding master edge
		bool bfound = true;
		for (int m=0; m<3; ++m)
		{
			if (fabs(E*Em[m]) > 0.9999)
			{
				FENodeSet& medge = masterEdges[m];

				int mref = closestNode(mesh, edge, rm);

				for (int i=0; i<(int)edge.size(); ++i)
				{
					assert(tag[edge[i]] < 0);
					if (tag[edge[i]] == -2)
					{
						vec3d ri = edge.Node(i)->m_r0;
						int k = closestNode(mesh, medge, ri);

						addLinearConstraint(*fem, edge[i], medge[k], edge[mref], refNode);
					}
				}
				
				bfound = true;
				break;
			}
		}
		assert(bfound);
	}

	return true;
}

/*
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
*/


int closestNode(FEMesh& mesh, const FENodeSet& set, const vec3d& r)
{
	int nmin = -1;
	double Dmin = 0.0;
	for (int i = 0; i<(int)set.size(); ++i)
	{
		vec3d& ri = mesh.Node(set[i]).m_r0;
		double D = (r - ri)*(r - ri);
		if ((D < Dmin) || (nmin == -1))
		{
			Dmin = D;
			nmin = i;
		}
	}
	return nmin;
}

// helper function for adding the linear constraints
void addLinearConstraint(FEModel& fem, int master, int slave, int nodeA, int nodeB)
{
	// get the linear constraint manager
	FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();

	// do one constraint for x, y, z
	for (int j = 0; j<3; ++j)
	{
		FELinearConstraint lc(&fem);
		lc.SetMasterDOF(j, master);

		lc.AddSlaveDof(j, slave, 1.0);
		lc.AddSlaveDof(j, nodeA, 1.0);
		lc.AddSlaveDof(j, nodeB, -1.0);

		LCM.AddLinearConstraint(lc);
	}
}
