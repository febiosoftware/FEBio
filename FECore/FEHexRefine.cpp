#include "stdafx.h"
#include "FEHexRefine.h"
#include "FEMesh.h"
#include "FEModel.h"
#include "FESolidDomain.h"
#include "FEMeshTopo.h"
#include "FEPrescribedDOF.h"
#include "FEElementList.h"
#include "FELinearConstraint.h"
#include "FELinearConstraintManager.h"
#include "log.h"

BEGIN_FECORE_CLASS(FEHexRefine, FEMeshAdaptor)
	ADD_PARAMETER(m_maxelem, "max_elems");
	ADD_PARAMETER(m_maxiter, "max_iter");

	ADD_PROPERTY(m_criterion, "criterion");
END_FECORE_CLASS();

FEHexRefine::FEHexRefine(FEModel* fem) : FERefineMesh(fem)
{
	m_maxelem = 0;
	m_maxiter = -1;
	m_criterion = nullptr;
}

bool FEHexRefine::Apply(int iteration)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	if ((m_maxiter >= 0) && (iteration >= m_maxiter))
	{
		feLog("\tMax iterations reached.\n");
		return true;
	}

	// see if we should do anything
	if ((m_maxelem > 0) && (mesh.Elements() >= m_maxelem))
	{
		feLog("\tElement limit reached.\n");
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
	// make sure we have a topo section
	if (m_topo == nullptr) return false;

	// refine the mesh
	if (RefineMesh(fem) == false)
	{
		feLog("\tNothing to do.\n");
		return false;
	}

	// update the BCs
	UpdateBCs(fem);

	// print some stats:
	FEMesh& mesh = fem.GetMesh();
	feLog("\tNodes : %d\n", mesh.Nodes());
	feLog("\tElements : %d\n", mesh.Elements());

	return true;
}

bool FEHexRefine::RefineMesh(FEModel& fem)
{
	FEMeshTopo& topo = *m_topo;

	FEMesh& mesh = fem.GetMesh();
	FEElementList allElems(mesh);

	const int NEL = mesh.Elements();
	const int NDOM = mesh.Domains();
	m_N0 = mesh.Nodes();
	m_NC = topo.Edges();
	int NF = topo.Faces();

	// Build the lists of items to split
	BuildSplitLists(fem);

	// make sure we have work to do
	if (m_splitElems == 0) return false;

	// we need to create a new node for each edge, face, and element that needs to be split
	int newNodes = m_splitEdges + m_splitFaces + m_splitElems;
	mesh.AddNodes(newNodes);

	// Next, the position and solution variables for all the nodes are updated.
	// Note that this has to be done before recreating the elements since 
	// the old elements are still needed to determine the new positions and solutions. 
	UpdateNewNodes(fem);

	// find the hanging nodes
	// and assign linear constraints to tie them down
	FindHangingNodes(fem);

	// Now, we can create new elements
	BuildNewDomains(fem);

	return true;
}

void FEHexRefine::BuildSplitLists(FEModel& fem)
{
	FEMeshTopo& topo = *m_topo;
	FEMesh& mesh = fem.GetMesh();

	// Get the elements that we need to refine
	int NEL = mesh.Elements();
	m_elemList.assign(NEL, -1);
	if (m_criterion)
	{
		vector<int> selection = m_criterion->GetElementList();
		for (int i = 0; i < selection.size(); ++i)
		{
			m_elemList[selection[i]] = 1;
		}
	}
	else
	{
		// just do'em all
		m_elemList.assign(NEL, 1);
	}

	// count how many elements to split
	int N1 = mesh.Nodes();
	m_splitElems = 0;
	for (int i = 0; i < m_elemList.size(); ++i) {
		if (m_elemList[i] == 1) {
			m_elemList[i] = N1++;
			m_splitElems++;
		}
	}

	// figure out which faces to refine
	int NF = topo.Faces();
	m_faceList.resize(NF, -1);
	FEElementFaceList& EFL = topo.m_EFL;
	for (int i = 0; i < m_elemList.size(); ++i)
	{
		if (m_elemList[i] != -1)
		{
			const std::vector<int>& elface = EFL.FaceList(i);
			for (int j = 0; j < elface.size(); ++j) m_faceList[elface[j]] = 1;
		}
	}

	// count how many faces to split
	m_splitFaces = 0;
	for (int i = 0; i < m_faceList.size(); ++i) {
		if (m_faceList[i] == 1) {
			m_faceList[i] = N1++;
			m_splitFaces++;
		}
	}

	// figure out which edges to refine
	m_edgeList.resize(m_NC, -1);
	FEFaceEdgeList& FEL = topo.m_FEL;
	for (int i = 0; i < NF; ++i)
	{
		if (m_faceList[i] != -1)
		{
			const std::vector<int>& faceEdge = FEL.EdgeList(i);
			for (int j = 0; j < faceEdge.size(); ++j) m_edgeList[faceEdge[j]] = 1;
		}
	}

	// count how many edges to split
	m_splitEdges = 0;
	for (int i = 0; i < m_NC; ++i) {
		if (m_edgeList[i] == 1) {
			m_edgeList[i] = N1++;
			m_splitEdges++;
		}
	}
}

void FEHexRefine::UpdateNewNodes(FEModel& fem)
{
	FEMeshTopo& topo = *m_topo;
	FEMesh& mesh = fem.GetMesh();

	// assign dofs to new nodes
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	m_NN = mesh.Nodes();
	for (int i = m_N0; i<m_NN; ++i)
	{
		FENode& node = mesh.Node(i);
		node.SetDOFS(MAX_DOFS);
	}

	// update the position of these new nodes
	for (int i = 0; i < topo.Edges(); ++i)
	{
		if (m_edgeList[i] != -1)
		{
			const FEEdgeList::EDGE& edge = topo.Edge(i);
			FENode& node = mesh.Node(m_edgeList[i]);
			vec3d r0 = mesh.Node(edge.node[0]).m_r0;
			vec3d r1 = mesh.Node(edge.node[1]).m_r0;
			node.m_r0 = (r0 + r1)*0.5;

			r0 = mesh.Node(edge.node[0]).m_rt;
			r1 = mesh.Node(edge.node[1]).m_rt;
			node.m_rt = (r0 + r1)*0.5;
		}
	}
	for (int i = 0; i < topo.Faces(); ++i)
	{
		if (m_faceList[i] != -1)
		{
			const FEFaceList::FACE& face = topo.Face(i);
			FENode& node = mesh.Node(m_faceList[i]);
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
	}
	int nelems = 0;
	FEElementList allElems(mesh);
	for (FEElementList::iterator it = allElems.begin(); it != allElems.end(); ++it, ++nelems)
	{
		if (m_elemList[nelems] != -1)
		{
			FESolidElement& el = dynamic_cast<FESolidElement&>(*it);

			int nn = el.Nodes();
			FENode& node = mesh.Node(m_elemList[nelems]);

			vec3d r0(0, 0, 0);
			for (int j = 0; j < nn; ++j) r0 += mesh.Node(el.m_node[j]).m_r0;
			r0 /= (double)nn;
			node.m_r0 = r0;

			vec3d rt(0, 0, 0);
			for (int j = 0; j < nn; ++j) rt += mesh.Node(el.m_node[j]).m_rt;
			rt /= (double)nn;
			node.m_rt = rt;
		}
	}

	// re-evaluate solution at nodes
	for (int i = 0; i < topo.Edges(); ++i)
	{
		if (m_edgeList[i] != -1)
		{
			const FEEdgeList::EDGE& edge = topo.Edge(i);
			FENode& node0 = mesh.Node(edge.node[0]);
			FENode& node1 = mesh.Node(edge.node[1]);

			FENode& node = mesh.Node(m_edgeList[i]);
			for (int j = 0; j < MAX_DOFS; ++j)
			{
				double v = (node0.get(j) + node1.get(j))*0.5;
				node.set(j, v);
			}
		}
	}
	for (int i = 0; i < topo.Faces(); ++i)
	{
		if (m_faceList[i] != -1)
		{
			const FEFaceList::FACE& face = topo.Face(i);
			FENode& node = mesh.Node(m_faceList[i]);
			int nn = face.ntype;
			for (int j = 0; j < MAX_DOFS; ++j)
			{
				double v = 0.0;
				for (int k = 0; k < nn; ++k) v += mesh.Node(face.node[k]).get(j);
				v /= (double)nn;
				node.set(j, v);
			}
		}
	}

	nelems = 0;
	for (FEElementList::iterator it = allElems.begin(); it != allElems.end(); ++it, ++nelems)
	{
		if (m_elemList[nelems] != -1)
		{
			FESolidElement& el = dynamic_cast<FESolidElement&>(*it);
			int nn = el.Nodes();
			FENode& node = mesh.Node(m_elemList[nelems]);
			for (int j = 0; j < MAX_DOFS; ++j)
			{
				double v = 0.0;
				for (int k = 0; k < nn; ++k) v += mesh.Node(el.m_node[k]).get(j);
				v /= (double)nn;
				node.set(j, v);
			}
		}
	}
}

void FEHexRefine::FindHangingNodes(FEModel& fem)
{
	FEMeshTopo& topo = *m_topo;
	FEMesh& mesh = fem.GetMesh();

	int hangingNodes = 0;

	FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();

	FEFaceEdgeList& FEL = topo.m_FEL;

	// we loop over non-split faces
	vector<int> tag(m_NC, 0);
	int NF = topo.Faces();
	for (int i = 0; i < NF; ++i)
	{
		if (m_faceList[i] == -1)
		{
			// This face is not split
			const std::vector<int>& fel = FEL.EdgeList(i);
			for (int j = 0; j < fel.size(); ++j)
			{
				if ((m_edgeList[fel[j]] >= 0) && (tag[fel[j]] == 0))
				{
					hangingNodes++;
					int nodeId = m_edgeList[fel[j]];

					// get the edge
					const FEEdgeList::EDGE& edge = topo.Edge(fel[j]);

					// setup a linear constraint for this node
					// TODO: This assumes mechanics! Need to generalize this to all active dofs!
					for (int k = 0; k < 3; ++k)
					{
						FELinearConstraint lc(&fem);
						lc.SetMasterDOF(k, nodeId);
						lc.AddSlaveDof(k, edge.node[0], 0.5);
						lc.AddSlaveDof(k, edge.node[1], 0.5);

						LCM.AddLinearConstraint(lc);
					}

					// set a tag to avoid double-counting
					tag[fel[j]] = 1;
				}
			}
		}
	}
	// we loop over the non-split elements
	int NEL = mesh.Elements();
	FEElementFaceList& EFL = topo.m_EFL;
	for (int i = 0; i<NEL; ++i)
	{
		if (m_elemList[i] == -1)
		{
			// This element is not split
			// If any of its faces are split, then the corresponding node
			// will be hanging
			const std::vector<int>& elface = EFL.FaceList(i);
			for (int j = 0; j < elface.size(); ++j)
			{
				if (m_faceList[elface[j]] >= 0)
				{
					hangingNodes++;

					int nodeId = m_faceList[elface[j]];

					// get the face
					const FEFaceList::FACE& face = topo.Face(elface[j]);

					// setup a linear constraint for this node
					// TODO: This assumes mechanics! Need to generalize this to all active dofs!
					for (int k = 0; k < 3; ++k)
					{
						FELinearConstraint lc(&fem);
						lc.SetMasterDOF(k, nodeId);
						lc.AddSlaveDof(k, face.node[0], 0.25);
						lc.AddSlaveDof(k, face.node[1], 0.25);
						lc.AddSlaveDof(k, face.node[2], 0.25);
						lc.AddSlaveDof(k, face.node[3], 0.25);

						LCM.AddLinearConstraint(lc);
					}

				}
			}
		}
	}

	if (hangingNodes > 0)
	{
		feLog("\tHanging nodes: %d\n", hangingNodes);
	}
}

void FEHexRefine::BuildNewDomains(FEModel& fem)
{
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

	FEMeshTopo& topo = *m_topo;

	FEMesh& mesh = fem.GetMesh();

	int nelems = 0;
	const int NDOM = mesh.Domains();
	for (int i = 0; i < NDOM; ++i)
	{
		// get the old domain
		FEDomain& oldDom = mesh.Domain(i);
		int NE0 = oldDom.Elements();

		// count how many elements to split in this domain
		int newElems = 0;
		for (int j = 0; j < NE0; ++j)
		{
			if (m_elemList[nelems + j] != -1) newElems++;
		}

		// make sure we have something to do
		if (newElems > 0)
		{
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
			oldDom.Create(8 * newElems + (NE0 - newElems), FE_HEX8G8);

			// set new element nodes
			int nel = 0;
			for (int j = 0; j < NE0; ++j, nelems++)
			{
				FEElement& el0 = newDom->ElementRef(j);

				if (m_elemList[nelems] != -1)
				{
					std::vector<int> ee = topo.m_EEL.EdgeList(nelems); assert(ee.size() == 12);
					std::vector<int> ef = topo.m_EFL.FaceList(nelems); assert(ef.size() == 6);

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
					ENL[8] = m_edgeList[ee[0]];
					ENL[9] = m_edgeList[ee[1]];
					ENL[10] = m_edgeList[ee[2]];
					ENL[11] = m_edgeList[ee[3]];
					ENL[12] = m_edgeList[ee[4]];
					ENL[13] = m_edgeList[ee[5]];
					ENL[14] = m_edgeList[ee[6]];
					ENL[15] = m_edgeList[ee[7]];
					ENL[16] = m_edgeList[ee[8]];
					ENL[17] = m_edgeList[ee[9]];
					ENL[18] = m_edgeList[ee[10]];
					ENL[19] = m_edgeList[ee[11]];
					ENL[20] = m_faceList[ef[0]];
					ENL[21] = m_faceList[ef[1]];
					ENL[22] = m_faceList[ef[2]];
					ENL[23] = m_faceList[ef[3]];
					ENL[24] = m_faceList[ef[4]];
					ENL[25] = m_faceList[ef[5]];
					ENL[26] = m_elemList[nelems];

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
				else
				{
					FEElement& el1 = oldDom.ElementRef(nel++);
					for (int k = 0; k < el0.Nodes(); ++k) el1.m_node[k] = el0.m_node[k];
				}
			}

			// we don't need this anymore
			delete newDom;
		}
	}

	// re-init domains
	for (int i = 0; i < NDOM; ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.CreateMaterialPointData();
		dom.Init();
		dom.Activate();
	}
}

void FEHexRefine::UpdateBCs(FEModel& fem)
{
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

	fem.GetLinearConstraintManager().Activate();
}
