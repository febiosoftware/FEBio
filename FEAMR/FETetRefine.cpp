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
#include "FETetRefine.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEMeshTopo.h>
#include <FECore/FEFixedBC.h>
#include <FECore/FESurface.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FETetRefine, FERefineMesh)
END_FECORE_CLASS();

FETetRefine::FETetRefine(FEModel* fem) : FERefineMesh(fem)
{
}

struct TRI
{
	int n[3];
};

bool FETetRefine::RefineMesh()
{
	FEModel& fem = *GetFEModel();

	if (m_topo == nullptr) return false;
	FEMeshTopo& topo = *m_topo;

	FEMesh& mesh = GetMesh();

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(0));
	int NEL = dom.Elements();

	// we need to create a new node for each edge
	int N0 = mesh.Nodes();
	int newNodes = topo.Edges();
	mesh.AddNodes(newNodes);
	int N1 = N0 + newNodes;

	// update the position of these new nodes
	int n = N0;
	for (int i = 0; i < topo.Edges(); ++i)
	{
		const FEEdgeList::EDGE& edge = topo.Edge(i);
		FENode& node = mesh.Node(n++);

		vec3d r0 = mesh.Node(edge.node[0]).m_r0;
		vec3d r1 = mesh.Node(edge.node[1]).m_r0;
		node.m_r0 = (r0 + r1)*0.5;

		r0 = mesh.Node(edge.node[0]).m_rt;
		r1 = mesh.Node(edge.node[1]).m_rt;
		node.m_rt = (r0 + r1)*0.5;
	}
	assert(n == N1);

	// assign dofs to new nodes
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	int NN = mesh.Nodes();
	for (int i = N0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		node.SetDOFS(MAX_DOFS);
	}

	// re-evaluate solution at nodes
	n = N0;
	for (int i = 0; i < topo.Edges(); ++i)
	{
		const FEEdgeList::EDGE& edge = topo.Edge(i);
		FENode& node0 = mesh.Node(edge.node[0]);
		FENode& node1 = mesh.Node(edge.node[1]);

		FENode& node = mesh.Node(n++);
		for (int j = 0; j < MAX_DOFS; ++j)
		{
			double v = (node0.get(j) + node1.get(j))*0.5;
			node.set(j, v);
		}
		node.UpdateValues();
	}
	assert(n == N1);

	const int LUT[8][4] = {
		{ 0, 4, 6, 7 },
		{ 4, 1, 5, 8 },
		{ 2, 6, 5, 9 },
		{ 7, 8, 9, 3 },
		{ 4, 8, 6, 7 },
		{ 4, 8, 5, 6 },
		{ 5, 6, 8, 9 },
		{ 6, 7, 8, 9 },
	};

	// now we recreate the domains
	const int NDOM = mesh.Domains();
	for (int i = 0; i < NDOM; ++i)
	{
		// get the old domain
		FEDomain& oldDom = mesh.Domain(i);
		int NE0 = oldDom.Elements();

		// create a copy of old domain (since we want to retain the old domain)
		FEDomain* newDom = fecore_new<FESolidDomain>(oldDom.GetTypeStr(), &fem);
		newDom->Create(NE0, FEElementLibrary::GetElementSpecFromType(FE_TET4G4));
		for (int j = 0; j < NE0; ++j)
		{
			FEElement& el0 = oldDom.ElementRef(j);
			FEElement& el1 = newDom->ElementRef(j);
			for (int k = 0; k < el0.Nodes(); ++k) el1.m_node[k] = el0.m_node[k];
		}

		// reallocate the old domain
		oldDom.Create(8 * NE0, FEElementLibrary::GetElementSpecFromType(FE_TET4G4));

		// set new element nodes
		int nel = 0;
		for (int j = 0; j < NE0; ++j)
		{
			FEElement& el0 = newDom->ElementRef(j);

			std::vector<int> ee = topo.ElementEdgeList(j); assert(ee.size() == 6);

			// build the look-up table
			int ENL[10] = { 0 };
			ENL[0] = el0.m_node[0];
			ENL[1] = el0.m_node[1];
			ENL[2] = el0.m_node[2];
			ENL[3] = el0.m_node[3];
			ENL[4] = N0 + ee[0];
			ENL[5] = N0 + ee[1];
			ENL[6] = N0 + ee[2];
			ENL[7] = N0 + ee[3];
			ENL[8] = N0 + ee[4];
			ENL[9] = N0 + ee[5];

			for (int k = 0; k < 8; ++k)
			{
				FEElement& el1 = oldDom.ElementRef(nel++);

				el1.m_node[0] = ENL[LUT[k][0]];
				el1.m_node[1] = ENL[LUT[k][1]];
				el1.m_node[2] = ENL[LUT[k][2]];
				el1.m_node[3] = ENL[LUT[k][3]];

				int a = 0;
			}
		}

		// we don't need this anymore
		delete newDom;
	}
	mesh.RebuildLUT();

	// re-init domains
	for (int i = 0; i < NDOM; ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.CreateMaterialPointData();
		dom.Reset();	// NOTE: we need to call this to actually call the Init function on the material points.
		dom.Init();
		dom.Activate();
	}

	const int FLUT[4][3] = {
		{ 0, 3, 5 },
		{ 1, 4, 3 },
		{ 2, 5, 4 },
		{ 3, 4, 5 },
	};

	// recreate element sets
	for (int i = 0; i < mesh.ElementSets(); ++i)
	{
		FEElementSet& eset = mesh.ElementSet(i);

		// get the domain list
		// NOTE: Don't get the reference, since then the same reference
		// is passed to Create below, which causes problems.
		FEDomainList domList = eset.GetDomainList();
		if (domList.IsEmpty()) { throw std::runtime_error("Error in FEMMGRemesh!"); }

		// recreate the element set from the domain list
		eset.Create(domList);
	}

	// recreate surfaces
	int faceMark = 1;
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		FESurface& surf = mesh.Surface(i);
		int NF0 = surf.Elements();

		vector<int> faceList = topo.FaceIndexList(surf);
		assert(faceList.size() == NF0);

		vector<TRI> tri(NF0);
		for (int j = 0; j < NF0; ++j)
		{
			FESurfaceElement& el = surf.Element(j);
			TRI t;
			t.n[0] = el.m_node[0];
			t.n[1] = el.m_node[1];
			t.n[2] = el.m_node[2];
			tri[j] = t;
		}

		int NF1 = NF0 * 4;
		surf.Create(NF1);
		int nf = 0;
		for (int j = 0; j < NF0; ++j)
		{
			TRI& t = tri[j];

			std::vector<int> ee = topo.FaceEdgeList(faceList[j]); assert(ee.size() == 3);

			// build the look-up table
			int FNL[6] = { 0 };
			FNL[0] = t.n[0];
			FNL[1] = t.n[1];
			FNL[2] = t.n[2];

			// the edges may not be ordered correctly, so we need to do a search here.
			for (int k = 0; k < 3; k++)
			{
				int a = t.n[k];
				int b = t.n[(k+1)%3];
				int e = -1;
				for (int l = 0; l < 3; ++l)
				{
					const FEEdgeList::EDGE& edge = topo.Edge(ee[l]);
					if ((edge.node[0] == a) && (edge.node[1] == b)) { e = l; break; }
					if ((edge.node[1] == a) && (edge.node[0] == b)) { e = l; break; }
				}
				assert(e >= 0);
				FNL[k + 3] = N0 + ee[e];
			}

			for (int k = 0; k < 4; ++k)
			{
				FESurfaceElement& fj = surf.Element(nf++);
				fj.SetType(FE_TRI3G3);
				fj.m_node[0] = FNL[FLUT[k][0]];
				fj.m_node[1] = FNL[FLUT[k][1]];
				fj.m_node[2] = FNL[FLUT[k][2]];
			}
		}

		surf.CreateMaterialPointData();
		surf.Init();

		// also update the facet set if the surface has one
		FEFacetSet* fset = surf.GetFacetSet();
		if (fset)
		{
			fset->Create(surf);
		}

		faceMark++;
	}

	// update node sets
	const int NSETS = mesh.NodeSets();
	for (int i=0; i<NSETS; ++i)
	{
		FENodeSet& nset = *mesh.NodeSet(i);
		vector<int> tag(mesh.Nodes(), 0);
		for (int j = 0; j < nset.Size(); ++j) tag[nset[j]] = 1;

		for (int j = 0; j < topo.Edges(); ++j)
		{
			const FEEdgeList::EDGE& edge = topo.Edge(j);

			if ((tag[edge.node[0]] == 1) && (tag[edge.node[1]] == 1))
			{
				nset.Add(N0 + j);
			}
		}
	}

	// re-activate the model
	UpdateModel();

	return true;
}
