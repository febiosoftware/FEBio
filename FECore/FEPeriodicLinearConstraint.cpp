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
#include "FEPeriodicLinearConstraint.h"
#include "FELinearConstraintManager.h"
#include "FELinearConstraint.h"
#include "FEModel.h"
#include "FESurface.h"

FEPeriodicLinearConstraint::FEPeriodicLinearConstraint(FEModel* fem) : FEStepComponent(fem), m_exclude(fem)
{
	m_refNode = -1;
}

FEPeriodicLinearConstraint::~FEPeriodicLinearConstraint()
{
}

void FEPeriodicLinearConstraint::AddNodeSetPair(const FENodeList& ms, const FENodeList& ss, bool push_back)
{
	NodeSetPair sp;
	sp.secondary = ms;
	sp.primary = ss;
	if (push_back) m_set.push_back(sp); else m_set.insert(m_set.begin(), sp);
}

// helper function for finding the closest node
int closestNode(FEMesh& mesh, const FENodeList& set, const vec3d& r);

// helper function for adding the linear constraints
void addLinearConstraint(FEModel& fem, int parent, int child, int nodeA, int nodeB);

// This function generates a linear constraint set based on the definition
// of three surface pairs.
bool FEPeriodicLinearConstraint::GenerateConstraints(FEModel* fem)
{
	// get the model's mesh
	FEMesh& mesh = GetMesh();

	// make sure there is a list of sets
	if (m_set.empty()) return true;

	// make sure we have three sets
	if (m_set.size() != 3) return false;

	// tag all nodes to identify what they are (edge, face, corner)
	int N = mesh.Nodes();
	vector<int> tag(N, 0);

	for (size_t i=0; i<m_set.size(); ++i)
	{
		FENodeList& ms = m_set[i].secondary;
		FENodeList& ss = m_set[i].primary;

		for (int j=0; j<(int)ms.Size(); ++j) tag[ms[j]]++;
		for (int j=0; j<(int)ss.Size(); ++j) tag[ss[j]]++;
	}

	// flip signs on primary
	for (size_t i = 0; i<m_set.size(); ++i)
	{
		FENodeList& ss = m_set[i].primary;
		for (int j = 0; j<(int)ss.Size(); ++j) 
		{
			int ntag = tag[ss[j]];
			if (ntag > 0) tag[ss[j]] = -ntag;
		}
	}

	// At this point, the following should hold
	// primary nodes: tag < 0, secondary nodes: tag > 0, interior nodes: tag = 0
	// only one secondary node should have a value of 3. We make this the reference node
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
		FENodeList& ms = m_set[n].secondary;
		FENodeList& ss = m_set[n].primary;

		// find the corresponding reference node on the primary surface
		int mref = closestNode(mesh, ss, rm);

		// make sure this is a corner node
		assert(tag[ss[mref]] == -3);

		// repeat for all primary nodes
		for (int i=0; i<(int)ss.Size(); ++i)
		{
			assert(tag[ss[i]] < 0);
			if (tag[ss[i]] == -1)
			{
				// get the primary node position
				vec3d& rs = ss.Node(i)->m_r0;

				// find the closest secondary node
				int m = closestNode(mesh, ms, rs);
				assert(tag[ms[m]] == 1);

				// add the linear constraints
				addLinearConstraint(*fem, ss[i], ms[m], ss[mref], refNode);
			}
		}
	}

	// extract all 12 edges
	vector<FENodeList> surf;
	for (int i=0; i<(int)m_set.size(); ++i)
	{
		surf.push_back(m_set[i].secondary);
		surf.push_back(m_set[i].primary);
	}

	vector<FENodeList> secondaryEdges;
	vector<FENodeList> primaryEdges;
	for (int i=0; i<surf.size(); ++i)
	{
		FENodeList& s0 = surf[i];
		for (int j=i+1; j<surf.size(); ++j)
		{
			FENodeList& s1 = surf[j];
			vector<int> tmp(N, 0);
			for (int k=0; k<s0.Size(); ++k) tmp[s0[k]]++;
			for (int k=0; k<s1.Size(); ++k) tmp[s1[k]]++;

			FENodeList edge(&mesh);
			for (int k=0; k<N; ++k)
			{
				if (tmp[k] == 2) edge.Add(k);
			}

			if (edge.Size() != 0)
			{
				// see if this is a secondary edge or not
				// we assume it's a secondary edge if it connects to the refnode
				bool bsecondary = false;
				for (int k=0; k<edge.Size(); ++k)
				{
					if (edge[k] == refNode)
					{
						bsecondary = true;
						break;
					}
				}

				if (bsecondary)
					secondaryEdges.push_back(edge);
				else
					primaryEdges.push_back(edge);
			}
		}
	}

	// since it is assumed the geometry is a cube, the following must hold
	assert(secondaryEdges.size() == 3);
	assert(primaryEdges.size() == 9);

	// find the secondary edge vectors
	vec3d Em[3];
	for (int i=0; i<3; ++i)
	{
		FENodeList& edge = secondaryEdges[i];

		// get the edge vector
		Em[i] = edge.Node(0)->m_r0 - edge.Node(1)->m_r0; assert(edge[0] != edge[1]); 
		Em[i].unit();
	}

	// setup linear constraints for edges
	for (int n = 0; n< primaryEdges.size(); ++n)
	{
		FENodeList& edge = primaryEdges[n];

		// get the edge vector
		vec3d E = edge.Node(0)->m_r0 - edge.Node(1)->m_r0; assert(edge[0] != edge[1]); E.unit();

		// find the corresponding secondary edge
		bool bfound = true;
		for (int m=0; m<3; ++m)
		{
			if (fabs(E*Em[m]) > 0.9999)
			{
				FENodeList& medge = secondaryEdges[m];

				int mref = closestNode(mesh, edge, rm);

				for (int i=0; i<(int)edge.Size(); ++i)
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

int closestNode(FEMesh& mesh, const FENodeList& set, const vec3d& r)
{
	int nmin = -1;
	double Dmin = 0.0;
	for (int i = 0; i<(int)set.Size(); ++i)
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
void addLinearConstraint(FEModel& fem, int parent, int child, int nodeA, int nodeB)
{
	// get the linear constraint manager
	FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();

	// do one constraint for x, y, z
	for (int j = 0; j<3; ++j)
	{
		FELinearConstraint* lc = fecore_alloc(FELinearConstraint, &fem);
		lc->SetParentDof(j, parent);

		lc->AddChildDof(j, child, 1.0);
		lc->AddChildDof(j, nodeA, 1.0);
		lc->AddChildDof(j, nodeB, -1.0);

		LCM.AddLinearConstraint(lc);
	}
}
