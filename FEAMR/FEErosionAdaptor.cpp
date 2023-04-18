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
#include "FEErosionAdaptor.h"
#include <FECore/FEMeshAdaptorCriterion.h>
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/FESurface.h>
#include <FECore/log.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FEElementList.h>
#include <FECore/FEMeshTopo.h>
#include <algorithm>
#include <stack>

BEGIN_FECORE_CLASS(FEErosionAdaptor, FEMeshAdaptor)
	ADD_PARAMETER(m_maxIters, "max_iters");
	ADD_PARAMETER(m_maxelem, "max_elems");
	ADD_PARAMETER(m_nsort, "sort");
	ADD_PARAMETER(m_bremoveIslands, "remove_islands");

	ADD_PROPERTY(m_criterion, "criterion");
END_FECORE_CLASS();

FEErosionAdaptor::FEErosionAdaptor(FEModel* fem) : FEMeshAdaptor(fem)
{
	m_maxIters = -1;
	m_maxelem = 0;
	m_nsort = 0;
	m_criterion = nullptr;
	m_bremoveIslands = false;
}

bool FEErosionAdaptor::Apply(int iteration)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	if ((m_maxIters >= 0) && (iteration >= m_maxIters))
	{
		feLog("\tMax iterations reached.");
		return false;
	}

	// Make sure there is a criterion
	if (m_criterion == nullptr) return false;

	// get the element selection
	FEMeshAdaptorSelection selection = m_criterion->GetElementSelection(GetElementSet());
	if (selection.empty())
	{
		feLog("\tNothing to do.\n");
		return false;
	}

	// see if we need to sort the data
	// (This is only necessary if the max_elem parameter is set)
	if ((m_maxelem > 0) && (m_nsort != 0))
	{
		if (m_nsort == 1)
		{
			// sort largest-to-smallest
			selection.Sort(FEMeshAdaptorSelection::SORT_DECREASING);
		}
		else if (m_nsort == 2)
		{
			// sort smallest-to-largest
			selection.Sort(FEMeshAdaptorSelection::SORT_INCREASING);
		}
	}

	// process the list
	int nsize = selection.size();
	if ((m_maxelem > 0) && (nsize > m_maxelem)) nsize = m_maxelem;
	for (int i = 0; i < nsize; ++i)
	{
		FEMeshAdaptorSelection::Item& it = selection[i];
		FEElement& el = *mesh.FindElementFromID(it.m_elementId);
		assert(el.isActive());
		el.setInactive();
	}
	feLog("\tDeactivated elements: %d\n", nsize);

	// remove any islands
	if (m_bremoveIslands) RemoveIslands();

	// if any nodes were orphaned, we need to deactivate them as well
	int NN = mesh.Nodes();
	vector<int> tag(NN, 0);
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j)
		{
			FEElement& el = dom.ElementRef(j);
			if (el.isActive())
			{
				int neln = el.Nodes();
				for (int n = 0; n < neln; ++n) tag[el.m_node[n]] = 1;
			}
		}
	}

	for (int i = 0; i < NN; ++i)
	{
		FENode& node = mesh.Node(i);
		if (tag[i] == 0)
		{
			node.SetFlags(FENode::EXCLUDE);
			int ndofs = node.dofs();
			for (int j = 0; j < ndofs; ++j)
				node.set_inactive(j);
		}
	}

	// any facets attached to a eroded element will be eroded as well. 
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		int erodedFaces = 0;
		FESurface& surf = mesh.Surface(i);
		for (int j = 0; j < surf.Elements(); ++j)
		{
			FESurfaceElement& face = surf.Element(j);
			FEElement* pe = face.m_elem[0]; assert(pe);
			if (pe && (pe->isActive() == false))
			{
				face.setInactive();
				erodedFaces++;
			}
		}
		if (erodedFaces != 0) surf.Init();
	}

	// remove any linear constraints of exclude nodes
	FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
	for (int j = 0; j < LCM.LinearConstraints();)
	{
		FELinearConstraint& lc = LCM.LinearConstraint(j);
		if (mesh.Node(lc.GetParentNode()).HasFlags(FENode::EXCLUDE))
		{
			LCM.RemoveLinearConstraint(j);
		}
		else ++j;
	}

	// also remove any linear constraints that have excluded child nodes
	for (int j = 0; j < LCM.LinearConstraints();)
	{
		FELinearConstraint& lc = LCM.LinearConstraint(j);

		bool del = false;
		int n = lc.Size();
		for (int k = 0; k < n; ++k)
		{
			if (mesh.Node(lc.GetChildDof(k).node).HasFlags(FENode::EXCLUDE))
			{
				del = true;
				break;
			}
		}
		if (del) LCM.RemoveLinearConstraint(j); else ++j;
	}

	// update model
	UpdateModel();

	return (nsize != 0);
}

void FEErosionAdaptor::RemoveIslands()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	FEMeshTopo topo;
	if (topo.Create(&mesh) == false)
	{
		feLogError("Failed removing islands.");
	}

	int NE = topo.Elements();
	vector<int> tag(NE, -1);

	// find an unprocessed element
	stack<int> s;
	int m = 0;	// island counter
	for (int n = 0; n < NE; ++n)
	{
		FEElement* el = topo.Element(n);
		if (el->isActive() && (tag[n] == -1))
		{
			// see if this is an island
			vector<int> island;

			tag[n] = m++;

			// push it on the stack
			s.push(n);
			while (s.empty() == false)
			{
				// pop the element
				int id = s.top(); s.pop();
				FEElement* el = topo.Element(id);
				island.push_back(id);

				// loop over all the neighbors
				vector<int> nbrList = topo.ElementNeighborIndexList(n);
				for (int i = 0; i < nbrList.size(); ++i)
				{
					FEElement* eli = topo.Element(nbrList[i]);
					if (eli && eli->isActive() && (tag[nbrList[i]] == -1))
					{
						tag[nbrList[i]] = m;
						s.push(nbrList[i]);
					}
				}
			}

			// Next, see if the island should be deactivated. 
			// It will be deactivated if all the nodes on the island are open
			bool isolated = true;
			for (int i = 0; i < island.size(); ++i)
			{
				FEElement* el = topo.Element(island[i]);

				int neln = el->Nodes();
				for (int j = 0; j < neln; ++j)
				{
					FENode& nj = mesh.Node(el->m_node[j]);

					// TODO: mechanics only!
					if ((nj.get_bc(0) != DOF_OPEN) ||
						(nj.get_bc(1) != DOF_OPEN) ||
						(nj.get_bc(2) != DOF_OPEN))
					{
						isolated = false;
						break;
					}
				}

				if (isolated == false) break;
			}

			if (isolated)
			{
				// island is isolated so deactivate all elements
				feLog("\tIsland of %d elements removed\n", island.size());
				for (int i = 0; i < island.size(); ++i)
				{
					FEElement* el = topo.Element(island[i]);
					el->setInactive();
				}
			}
		}
	}
}
