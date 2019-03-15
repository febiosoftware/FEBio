#include "stdafx.h"
#include "FEHexRefine.h"
#include "FEMesh.h"
#include "FEModel.h"
#include "FESolidDomain.h"
#include "FEMeshTopo.h"
#include "FEPrescribedDOF.h"
#include "log.h"

BEGIN_FECORE_CLASS(FEHexRefine, FEMeshAdaptor)
	ADD_PARAMETER(m_maxelem, "max_elem");
END_FECORE_CLASS();

FEHexRefine::FEHexRefine(FEModel* fem) : FERefineMesh(fem)
{
	m_maxelem = 0;
}

bool FEHexRefine::Apply(int iteration)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// see if we should do anything
	if ((m_maxelem > 0) && (mesh.Elements() >= m_maxelem))
	{
		feLog("Element limit reached.\n");
		return true;
	}

	// make sure this is a hex mesh
	// Note that we return true on error to indicate that 
	// the mesh did not change
	// TODO: Maybe trhow an exception instead of just returning on error
	if (mesh.IsType(ET_HEX8) == false)
	{
		feLogError("Cannot apply hex refinement: Mesh is not a HEX8 mesh.");
		return true;
	}

	// build the mesh-topo
	// Note that we return true on error to indicate that 
	// the mesh did not change
	// TODO: Maybe trhow an exception instead of just returning on error
	if (BuildMeshTopo(fem) == false)
	{
		feLogError("Cannot apply hex refinement: Error building topo structure.");
		return true;
	}

	// do the mesh refinement
	return !DoHexRefinement(fem);
}

bool FEHexRefine::DoHexRefinement(FEModel& fem)
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
		{ 0,  8, 24, 11, 16, 20, 26, 23 },
		{ 8,  1,  9, 24, 20, 17, 21, 26 },
		{ 11, 24, 10,  3, 23, 26, 22, 19 },
		{ 24,  9,  2, 10, 26, 21, 18, 22 },
		{ 16, 20, 26, 23,  4, 12, 25, 15 },
		{ 20, 17, 21, 26, 12,  5, 13, 25 },
		{ 23, 26, 22, 19, 15, 25, 14,  7 },
		{ 26, 21, 18, 22, 25, 13,  6, 14 },
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
			ENL[0] = el0.m_node[0];
			ENL[1] = el0.m_node[1];
			ENL[2] = el0.m_node[2];
			ENL[3] = el0.m_node[3];
			ENL[4] = el0.m_node[4];
			ENL[5] = el0.m_node[5];
			ENL[6] = el0.m_node[6];
			ENL[7] = el0.m_node[7];
			ENL[8] = m_N0 + ee[0];
			ENL[9] = m_N0 + ee[1];
			ENL[10] = m_N0 + ee[2];
			ENL[11] = m_N0 + ee[3];
			ENL[12] = m_N0 + ee[4];
			ENL[13] = m_N0 + ee[5];
			ENL[14] = m_N0 + ee[6];
			ENL[15] = m_N0 + ee[7];
			ENL[16] = m_N0 + ee[8];
			ENL[17] = m_N0 + ee[9];
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