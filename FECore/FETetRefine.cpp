#include "stdafx.h"
#include "FETetRefine.h"
#include "FESolidDomain.h"
#include "FEMeshTopo.h"
#include "FEFixedBC.h"

FETetRefine::FETetRefine(FEModel* fem) : FERefineMesh(fem)
{

}

bool FETetRefine::Apply(int iteration)
{
	FEModel& fem = *GetFEModel();

	// build the mesh-topo
	if (BuildMeshTopo(fem) == false) return false;

	// do the mesh refinement
	if (DoTetRefinement(fem) == false) return false;

	// all done
	return true;
}

bool FETetRefine::DoTetRefinement(FEModel& fem)
{
	if (m_topo == nullptr) return false;
	FEMeshTopo& topo = *m_topo;

	FEMesh& mesh = fem.GetMesh();

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(0));
	int NEL = dom.Elements();

	// we need to create a new node for each edge
	int N0 = mesh.Nodes();
	int newNodes = topo.m_edgeList.Edges();
	mesh.AddNodes(newNodes);
	int N1 = N0 + newNodes;

	// update the position of these new nodes
	int n = N0;
	for (int i = 0; i < topo.m_edgeList.Edges(); ++i)
	{
		const FEEdgeList::EDGE& edge = topo.m_edgeList.Edge(i);
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
	m_NN = mesh.Nodes();
	for (int i = N0; i<m_NN; ++i)
	{
		FENode& node = mesh.Node(i);
		node.SetDOFS(MAX_DOFS);
	}

	// re-evaluate solution at nodes
	n = N0;
	for (int i = 0; i < topo.m_edgeList.Edges(); ++i)
	{
		const FEEdgeList::EDGE& edge = topo.m_edgeList.Edge(i);
		FENode& node0 = mesh.Node(edge.node[0]);
		FENode& node1 = mesh.Node(edge.node[1]);

		FENode& node = mesh.Node(n++);
		for (int j = 0; j < MAX_DOFS; ++j)
		{
			double v = (node0.get(j) + node1.get(j))*0.5;
			node.set(j, v);
		}
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
		FEDomain* newDom = fecore_new<FEDomain>(oldDom.GetTypeStr(), &fem);
		newDom->Create(NE0, FE_TET4G4);
		for (int j = 0; j < NE0; ++j)
		{
			FEElement& el0 = oldDom.ElementRef(j);
			FEElement& el1 = newDom->ElementRef(j);
			for (int k = 0; k < el0.Nodes(); ++k) el1.m_node[k] = el0.m_node[k];
		}

		// reallocate the old domain
		oldDom.Create(8 * NE0, FE_TET4G4);

		// set new element nodes
		int nel = 0;
		for (int j = 0; j < NE0; ++j)
		{
			FEElement& el0 = newDom->ElementRef(j);

			std::vector<int> ee = topo.m_EEL.EdgeList(j); assert(ee.size() == 6);

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

	// re-init domains
	for (int i = 0; i < NDOM; ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.CreateMaterialPointData();
		dom.Init();
	}

	// translate BCs
	vector<int> tag(m_NN, 0);
	for (int i = 0; i < fem.FixedBCs(); ++i)
	{
		FEFixedBC& bc = *fem.FixedBC(i);

		vector<int> nodeList = bc.GetNodeList();
		for (int j = 0; j < (int)nodeList.size(); ++j) tag[nodeList[j]] = 1;

		for (int j = 0; j < topo.m_edgeList.Edges(); ++j)
		{
			const FEEdgeList::EDGE& edge = topo.m_edgeList[j];

			if ((tag[edge.node[0]] == 1) && (tag[edge.node[1]] == 1))
			{
				nodeList.push_back(N0 + j);
			}
		}

		// set the node list
		bc.SetNodeList(nodeList);

		// re-activate the bc
		if (bc.IsActive()) bc.Activate();
	}

	return true;
}
