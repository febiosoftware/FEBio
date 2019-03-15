#include "stdafx.h"
#include "FERefineMesh.h"
#include "FEModel.h"
#include "FESolidDomain.h"
#include "FEEdgeList.h"
#include "FEElementList.h"
#include "FEFaceList.h"
#include "FEFixedBC.h"
#include "FEPrescribedDOF.h"
#include "FEMeshTopo.h"
#include "log.h"

BEGIN_FECORE_CLASS(FERefineMesh, FEMeshAdaptor)
	ADD_PARAMETER(m_doMeshReset, "reset_mesh");
	ADD_PARAMETER(m_maxelem, "max_elem");
END_FECORE_CLASS();

FERefineMesh::FERefineMesh(FEModel* fem) : FEMeshAdaptor(fem), m_topo(nullptr)
{
	m_doMeshReset = false;
	m_maxelem = 0;
}

void FERefineMesh::DoMeshReset(bool b)
{
	m_doMeshReset = b;
}

bool FERefineMesh::BuildMeshTopo(FEModel& fem)
{
	if (m_topo) { delete m_topo; m_topo = nullptr; }
	m_topo = new FEMeshTopo;
	return m_topo->Create(&fem.GetMesh());
}

bool FERefineMesh::Apply()
{
	FEModel& fem = *GetFEModel();

	// build the mesh-topo
	if (BuildMeshTopo(fem) == false) return false;

	// do the mesh refinement
	if (DoMeshRefinement(fem) == false) return false;

	// Reset the model if requested
	if (m_doMeshReset)
	{
		// we only reset the mesh and then re-activate the boundary conditions
		fem.GetMesh().Reset();
		fem.Activate();
	}

	// all done
	return true;
}

// Here, the return code indicates whether the time step should be rerun or not
// So, return yes if the mesh was not modified, false otherwise
bool FERefineMesh::DoMeshRefinement(FEModel& fem)
{
	// see if this model is a hex or tet (or neither)
	FEMesh& mesh = fem.GetMesh();

	// see if we should do anything
	if ((m_maxelem > 0) && (mesh.Elements() >= m_maxelem))
	{
		feLog("Element limit reached.\n");
		return true;
	}

	if (mesh.IsType(ET_TET4))
	{
		return !DoTetRefinement(fem);
	}
	else if (mesh.IsType(ET_HEX8))
	{
		return !DoHexRefinement(fem);
	}

	return true;
}

bool FERefineMesh::DoTetRefinement(FEModel& fem)
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
		{0, 4, 6, 7},
		{4, 1, 5, 8},
		{2, 6, 5, 9},
		{7, 8, 9, 3},
		{4, 8, 6, 7},
		{4, 8, 5, 6},
		{5, 6, 8, 9},
		{6, 7, 8, 9},
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
	m_tag.resize(m_NN, 0);
	for (int i = 0; i < fem.FixedBCs(); ++i)
	{
		FEFixedBC& bc = *fem.FixedBC(i);

		vector<int> nodeList = bc.GetNodeList();
		for (int j = 0; j < (int)nodeList.size(); ++j) m_tag[nodeList[j]] = 1;

		for (int j = 0; j < topo.m_edgeList.Edges(); ++j)
		{
			const FEEdgeList::EDGE& edge = topo.m_edgeList[j];

			if ((m_tag[edge.node[0]] == 1) && (m_tag[edge.node[1]] == 1))
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

bool FERefineMesh::DoHexRefinement(FEModel& fem)
{
	if (m_topo == nullptr) return false;
	FEMeshTopo& topo = *m_topo;

	FEMesh& mesh = fem.GetMesh();

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(0));
	int NEL = dom.Elements();

	// we need to create a new node for each edge, face, and element
	m_N0 = mesh.Nodes();
	m_NC = topo.m_edgeList.Edges();
	int NF = topo.m_faceList.Faces();
	int newNodes = m_NC + NF + NEL;
	mesh.AddNodes(newNodes);
	int N1 = m_N0 + newNodes;

	// assign dofs to new nodes
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	m_NN = mesh.Nodes();
	for (int i = m_N0; i<m_NN; ++i)
	{
		FENode& node = mesh.Node(i);
		node.SetDOFS(MAX_DOFS);
	}

	// Now, I update the position and solution variables for all the nodes
	// Note that I do this before recreating the elements since I still need
	// the old elements to determine the new positions and solutions. 

	// update the position of these new nodes
	int n = m_N0;
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
	for (int i = 0; i < topo.m_faceList.Faces(); ++i)
	{
		const FEFaceList::FACE& face = topo.m_faceList.Face(i);
		FENode& node = mesh.Node(n++);
		int nn = face.ntype;

		vec3d r0(0, 0, 0);
		for (int j = 0; j < nn; ++j) r0 += mesh.Node(face.node[j]).m_r0;
		r0 /= (double)nn;
		node.m_r0 = r0;

		vec3d rt(0, 0, 0);
		for (int j = 0; j < nn; ++j) rt += mesh.Node(face.node[j]).m_rt;
		rt /= (double)nn;
		node.m_rt = rt;
	}
	for (int i = 0; i < NEL; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nn = el.Nodes();
		FENode& node = mesh.Node(n++);

		vec3d r0(0, 0, 0);
		for (int j = 0; j < nn; ++j) r0 += mesh.Node(el.m_node[j]).m_r0;
		r0 /= (double)nn;
		node.m_r0 = r0;

		vec3d rt(0, 0, 0);
		for (int j = 0; j < nn; ++j) rt += mesh.Node(el.m_node[j]).m_rt;
		rt /= (double)nn;
		node.m_rt = rt;
	}
	assert(n == N1);

	// re-evaluate solution at nodes
	n = m_N0;
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
	for (int i = 0; i < topo.m_faceList.Faces(); ++i)
	{
		const FEFaceList::FACE& face = topo.m_faceList.Face(i);
		FENode& node = mesh.Node(n++);
		int nn = face.ntype;
		for (int j = 0; j < MAX_DOFS; ++j)
		{
			double v = 0.0;
			for (int k = 0; k < nn; ++k) v += mesh.Node(face.node[k]).get(j);
			v /= (double)nn;
			node.set(j, v);
		}
	}
	for (int i = 0; i < NEL; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nn = el.Nodes();
		FENode& node = mesh.Node(n++);
		for (int j = 0; j < MAX_DOFS; ++j)
		{
			double v = 0.0;
			for (int k = 0; k < nn; ++k) v += mesh.Node(el.m_node[k]).get(j);
			v /= (double)nn;
			node.set(j, v);
		}
	}
	assert(n == N1);

	// Now, we can create new elements
	const int LUT[8][8] = {
		{  0,  8, 24, 11, 16, 20, 26, 23},
		{  8,  1,  9, 24, 20, 17, 21, 26},
		{ 11, 24, 10,  3, 23, 26, 22, 19},
		{ 24,  9,  2, 10, 26, 21, 18, 22},
		{ 16, 20, 26, 23,  4, 12, 25, 15},
		{ 20, 17, 21, 26, 12,  5, 13, 25},
		{ 23, 26, 22, 19, 15, 25, 14,  7},
		{ 26, 21, 18, 22, 25, 13,  6, 14},
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
		newDom->Create(NE0, FE_HEX8G8);
		for (int j = 0; j < NE0; ++j)
		{
			FEElement& el0 = oldDom.ElementRef(j);
			FEElement& el1 = newDom->ElementRef(j);
			for (int k = 0; k < el0.Nodes(); ++k) el1.m_node[k] = el0.m_node[k];
		}

		// reallocate the old domain
		oldDom.Create(8 * NE0, FE_HEX8G8);

		// set new element nodes
		int nel = 0;
		for (int j = 0; j < NE0; ++j)
		{
			FEElement& el0 = newDom->ElementRef(j);

			std::vector<int> ee = topo.m_EEL.EdgeList(j); assert(ee.size() == 12);
			std::vector<int> ef = topo.m_EFL.FaceList(j); assert(ef.size() == 6);

			// build the look-up table
			int ENL[27] = { 0 };
			ENL[ 0] = el0.m_node[0];
			ENL[ 1] = el0.m_node[1];
			ENL[ 2] = el0.m_node[2];
			ENL[ 3] = el0.m_node[3];
			ENL[ 4] = el0.m_node[4];
			ENL[ 5] = el0.m_node[5];
			ENL[ 6] = el0.m_node[6];
			ENL[ 7] = el0.m_node[7];
			ENL[ 8] = m_N0 + ee[ 0];
			ENL[ 9] = m_N0 + ee[ 1];
			ENL[10] = m_N0 + ee[ 2];
			ENL[11] = m_N0 + ee[ 3];
			ENL[12] = m_N0 + ee[ 4];
			ENL[13] = m_N0 + ee[ 5];
			ENL[14] = m_N0 + ee[ 6];
			ENL[15] = m_N0 + ee[ 7];
			ENL[16] = m_N0 + ee[ 8];
			ENL[17] = m_N0 + ee[ 9];
			ENL[18] = m_N0 + ee[10];
			ENL[19] = m_N0 + ee[11];
			ENL[20] = m_N0 + m_NC + ef[0];
			ENL[21] = m_N0 + m_NC + ef[1];
			ENL[22] = m_N0 + m_NC + ef[2];
			ENL[23] = m_N0 + m_NC + ef[3];
			ENL[24] = m_N0 + m_NC + ef[4];
			ENL[25] = m_N0 + m_NC + ef[5];
			ENL[26] = m_N0 + m_NC + NF + j;

			for (int k = 0; k < 8; ++k)
			{
				FEElement& el1 = oldDom.ElementRef(nel++);

				el1.m_node[0] = ENL[LUT[k][0]];
				el1.m_node[1] = ENL[LUT[k][1]];
				el1.m_node[2] = ENL[LUT[k][2]];
				el1.m_node[3] = ENL[LUT[k][3]];
				el1.m_node[4] = ENL[LUT[k][4]];
				el1.m_node[5] = ENL[LUT[k][5]];
				el1.m_node[6] = ENL[LUT[k][6]];
				el1.m_node[7] = ENL[LUT[k][7]];
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
		dom.Activate();
	}

	// translate BCs
	for (int i = 0; i < fem.FixedBCs(); ++i)
	{
		FEFixedBC& bc = *fem.FixedBC(i);
		UpdateFixedBC(bc);
	}

	for (int i = 0; i < fem.PrescribedBCs(); ++i)
	{
		FEPrescribedDOF& bc = dynamic_cast<FEPrescribedDOF&>(*fem.PrescribedBC(i));
		UpdatePrescribedBC(bc);
	}


	// print some stats:
	feLog("Nodes : %d\n", mesh.Nodes());
	feLog("Elements : %d\n", mesh.Elements());

	return true;
}

void FERefineMesh::UpdateFixedBC(FEFixedBC& bc)
{
	vector<int> nodeList = bc.GetNodeList();

	FEMeshTopo& topo = *m_topo;
	vector<int> tag(m_NN, 0);

	for (int j = 0; j < (int)nodeList.size(); ++j) tag[nodeList[j]] = 1;

	for (int j = 0; j < topo.m_edgeList.Edges(); ++j)
	{
		const FEEdgeList::EDGE& edge = topo.m_edgeList[j];

		if ((tag[edge.node[0]] == 1) && (tag[edge.node[1]] == 1))
		{
			nodeList.push_back(m_N0 + j);
		}
	}

	for (int j = 0; j < topo.m_faceList.Faces(); ++j)
	{
		const FEFaceList::FACE& face = topo.m_faceList.Face(j);

		assert(face.ntype == 4);
		if ((tag[face.node[0]] == 1) &&
			(tag[face.node[1]] == 1) &&
			(tag[face.node[2]] == 1) &&
			(tag[face.node[3]] == 1))
		{
			nodeList.push_back(m_N0 + m_NC + j);
		}
	}

	// set the node list
	bc.SetNodeList(nodeList);

	// re-activate the bc
	if (bc.IsActive()) bc.Activate();
}

void FERefineMesh::UpdatePrescribedBC(FEPrescribedDOF& bc)
{
	int items = bc.Items();

	FEMeshTopo& topo = *m_topo;
	vector<int> tag(m_NN, -1);

	for (int j = 0; j < items; ++j) tag[bc.GetItem(j).nid] = j;

	for (int j = 0; j < topo.m_edgeList.Edges(); ++j)
	{
		const FEEdgeList::EDGE& edge = topo.m_edgeList[j];

		if ((tag[edge.node[0]] >= 0) && (tag[edge.node[1]] >= 0))
		{
			double a0 = bc.GetItem(tag[edge.node[0]]).ref;
			double a1 = bc.GetItem(tag[edge.node[1]]).ref;
			bc.AddNode(m_N0 + j, (a0 + a1)*0.5);
		}
	}

	for (int j = 0; j < topo.m_faceList.Faces(); ++j)
	{
		const FEFaceList::FACE& face = topo.m_faceList.Face(j);

		assert(face.ntype == 4);
		if ((tag[face.node[0]] >= 0) &&
			(tag[face.node[1]] >= 0) &&
			(tag[face.node[2]] >= 0) &&
			(tag[face.node[3]] >= 0))
		{
			double a0 = bc.GetItem(tag[face.node[0]]).ref;
			double a1 = bc.GetItem(tag[face.node[1]]).ref;
			double a2 = bc.GetItem(tag[face.node[2]]).ref;
			double a3 = bc.GetItem(tag[face.node[3]]).ref;
			bc.AddNode(m_N0 + m_NC + j, (a0 + a1 + a2 + a3)*0.25);
		}
	}

	// re-activate the bc
	if (bc.IsActive()) bc.Activate();
}
