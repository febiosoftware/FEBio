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
#include <set>

BEGIN_FECORE_CLASS(FEErosionAdaptor, FEMeshAdaptor)
	ADD_PARAMETER(m_maxIters, "max_iters");
	ADD_PARAMETER(m_maxelem, "max_elems");
	ADD_PARAMETER(m_nsort, "sort");
	ADD_PARAMETER(m_bremoveIslands, "remove_islands");
	ADD_PARAMETER(m_erodeSurfaces, "erode_surfaces")->setEnums("no\0yes\0grow\0reconstruct\0");

	ADD_PROPERTY(m_criterion, "criterion");
END_FECORE_CLASS();

FEErosionAdaptor::FEErosionAdaptor(FEModel* fem) : FEMeshAdaptor(fem)
{
	m_maxIters = -1;
	m_maxelem = 0;
	m_nsort = 0;
	m_erodeSurfaces = SurfaceErodeOption::ERODE;

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
	int deactivatedElements = 0;
	for (int i = 0; i < nsize; ++i)
	{
		FEMeshAdaptorSelection::Item& it = selection[i];
		FEElement& el = *mesh.FindElementFromID(it.m_elementId);
		if (el.isActive())
		{
			el.setInactive();
			deactivatedElements++;
		}
	}
	feLog("\tDeactivated elements: %d\n", deactivatedElements);
	if (deactivatedElements == 0) return false;

	FEMeshTopo topo;
	if (topo.Create(&mesh) == false)
	{
		feLogError("Failed building mesh topo.");
		return false;
	}

	// remove any islands
	if (m_bremoveIslands) RemoveIslands(topo);

	// if any nodes were orphaned, we need to deactivate them as well
	DeactivateOrphanedNodes();

	// any facets attached to a eroded element will be eroded as well.
	switch (m_erodeSurfaces)
	{
	case SurfaceErodeOption::DONT_ERODE: break;// don't do anything
	case SurfaceErodeOption::ERODE: ErodeSurfaces(); break;
	case SurfaceErodeOption::GROW : GrowErodedSurfaces(topo); break;
	case SurfaceErodeOption::RECONSTRUCT: ReconstructSurfaces(topo); break;
	}

	// remove any linear constraints of excluded nodes
	UpdateLinearConstraints();

	// update model
	UpdateModel();

	return (nsize != 0);
}

void FEErosionAdaptor::RemoveIslands(FEMeshTopo& topo)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

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
				vector<int> nbrList = topo.ElementNeighborIndexList(id);
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

void FEErosionAdaptor::DeactivateOrphanedNodes()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

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
}

void FEErosionAdaptor::ErodeSurfaces()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		int erodedFaces = 0;
		FESurface& surf = mesh.Surface(i);
		for (int j = 0; j < surf.Elements(); ++j)
		{
			FESurfaceElement& face = surf.Element(j);
			FEElement* pe = face.m_elem[0].pe; assert(pe);
			if (pe && (pe->isActive() == false))
			{
				face.setInactive();
				erodedFaces++;
			}
		}
		if (erodedFaces != 0) surf.Init();
	}
}

void FEErosionAdaptor::GrowErodedSurfaces(FEMeshTopo& topo)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		FESurface& surf = mesh.Surface(i);

		// do a search to see if we need to modify this surface
		bool doErode = false;
		for (int j = 0; j < surf.Elements(); ++j)
		{
			FESurfaceElement& face = surf.Element(j);
			if (face.isActive())
			{
				FEElement* pe = face.m_elem[0].pe; assert(pe);
				if (pe && (pe->isActive() == false))
				{
					doErode = true;
					break;
				}
			}
		}

		if (doErode)
		{
			std::vector<FEFaceList::FACE> newFaces;
			std::set<FEElement*> processedElements;
			newFaces.reserve(surf.Elements());
			int erodedFaces = 0;
			for (int j = 0; j < surf.Elements(); ++j)
			{
				FESurfaceElement& face = surf.Element(j);
				if (face.isActive())
				{
					FEElement* pe = face.m_elem[0].pe; assert(pe);
					if (pe && (pe->isActive() == false))
					{
						// make sure we didn't process this element yet
						if (processedElements.find(pe) == processedElements.end())
						{
							int id = topo.GetElementIndexFromID(pe->GetID());
							std::vector<FEElement*> nbrList = topo.ElementNeighborList(id);
							for (int k = 0; k < nbrList.size(); ++k)
							{
								FEElement* pk = nbrList[k];
								if ((pk != nullptr) && pk->isActive())
								{
									int node[FEElement::MAX_NODES] = { 0 };
									int nf = pe->GetFace(k, node);

									// Note that we invert the element to
									// make sure the new face is outward pointing
									FEFaceList::FACE newFace;
									newFace.ntype = nf;
									for (int l = 0; l < nf; ++l)
										newFace.node[l] = node[nf - 1 - l];

									if (!newFace.IsEqual(face.m_node.data()))
										newFaces.push_back(newFace);
								}
							}

							processedElements.insert(pe);
						}

						face.setInactive();
						erodedFaces++;
					}
					else
					{
						FEFaceList::FACE newFace;
						newFace.ntype = face.Nodes();
						for (int k = 0; k < face.Nodes(); ++k)
							newFace.node[k] = face.m_node[k];
						newFaces.push_back(newFace);
					}
				}
			}
			feLog("eroded faces: %d", erodedFaces);
			if (erodedFaces != 0)
			{
				int faces = newFaces.size();
				surf.Create(faces);
				for (int j = 0; j < faces; ++j)
				{
					FEFaceList::FACE& f = newFaces[j];
					FESurfaceElement& face = surf.Element(j);
					face.setActive();

					if (f.ntype == 3)
						face.SetType(FE_TRI3G3);
					else if (f.ntype == 4)
						face.SetType(FE_QUAD4G4);

					for (int k = 0; k < f.ntype; ++k)
						face.m_node[k] = f.node[k];
				}

				surf.CreateMaterialPointData();
				surf.Init();

				// also update the facet set if the surface has one
				FEFacetSet* fset = surf.GetFacetSet();
				if (fset)
				{
					fset->Create(surf);
				}
			}
		}
	}
}

void FEErosionAdaptor::UpdateLinearConstraints()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

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
}

void FEErosionAdaptor::ReconstructSurfaces(FEMeshTopo& topo)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// We assume that the surfaces are definded through part lists. 
	// We don't actually store which surfaces are generated from part lists,
	// but they should have the same name. 
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		FESurface& surf = mesh.Surface(i);

		// see if a part list exists with this name
		FEDomainList* domList = mesh.FindDomainList(surf.GetName());
		if (domList)
		{
			// build the outer surface
			FEFacetSet* facetSet = mesh.DomainBoundary(*domList);

			FEFacetSet* oldFacetSet = surf.GetFacetSet();
			if (oldFacetSet == nullptr)
			{
				surf.Create(*facetSet);
			}
			else
			{
				oldFacetSet->Clear();
				oldFacetSet->Add(facetSet);
				surf.Create(*oldFacetSet);
			}

			surf.CreateMaterialPointData();
			surf.Init();
		}
	}
}
