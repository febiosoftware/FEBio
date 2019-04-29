/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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

BEGIN_FECORE_CLASS(FEHexRefine, FERefineMesh)
	ADD_PARAMETER(m_maxelem, "max_elems");
	ADD_PARAMETER(m_maxiter, "max_iter");
	ADD_PARAMETER(m_elemRefine, "max_elem_refine");
	ADD_PROPERTY(m_criterion, "criterion");
END_FECORE_CLASS();

FEHexRefine::FEHexRefine(FEModel* fem) : FERefineMesh(fem)
{
	m_maxelem = 0;
	m_elemRefine = 0;
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
	if (BuildMeshTopo() == false)
	{
		feLogError("Cannot apply hex refinement: Error building topo structure.");
		return true;
	}

	// refine the mesh
	if (RefineMesh(fem) == false)
	{
		feLog("\tNothing to do.\n");
		return true;
	}

	// update the BCs
	// TODO: Should I deactivate all BCs prior to doing the mesh refinement?
	UpdateBCs();

	return false;
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

	// Next, the position and solution variables for all the nodes are updated.
	// Note that this has to be done before recreating the elements since 
	// the old elements are still needed to determine the new positions and solutions. 
	UpdateNewNodes(fem);

	// find the hanging nodes
	// and assign linear constraints to tie them down
	FindHangingNodes(fem);

	// Now, we can create new elements
	BuildNewDomains(fem);

	// print some stats:
	feLog("\tNew mesh stats:\n");
	feLog("\t  Nodes .......... : %d\n", mesh.Nodes());
	feLog("\t  Elements ....... : %d\n", mesh.Elements());

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

	// We cannot split elements that have a hanging nodes
	// so remove those elements from the list
	int nrejected = 0;
	for (int i = 0; i < NEL; ++i)
	{
		if (m_elemList[i] != -1)
		{
			FEElement& el = *topo.Element(i);
			int nel = el.Nodes();
			for (int j = 0; j < nel; ++j)
			{
				if (mesh.Node(el.m_node[j]).HasFlags(FENode::HANGING))
				{
					// sorry, can't split this element
					m_elemList[i] = -1;
					nrejected++;
					break;
				}
			}
		}
	}

	if (nrejected > 0)
	{
		feLog("\tElements rejected: %d\n", nrejected);
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
	if (m_splitElems == 0) return;

	// make sure we don't exceed the max elements per refinement step
	if ((m_elemRefine > 0) && (m_splitElems > m_elemRefine))
	{
		N1 = mesh.Nodes();
		m_splitElems = 0;
		for (int i = 0; i < m_elemList.size(); ++i) {
			if (m_elemList[i] >= 0) {
				if (m_splitElems >= m_elemRefine)
				{
					m_elemList[i] = -1;
				}
				else
				{
					m_splitElems++;
					m_elemList[i] = N1++;
				}
			}
		}
		assert(m_splitElems == m_elemRefine);
	}

	// figure out which faces to refine
	int NF = topo.Faces();
	m_faceList.assign(NF, -1);
	for (int i = 0; i < m_elemList.size(); ++i)
	{
		if (m_elemList[i] != -1)
		{
			const std::vector<int>& elface = topo.ElementFaceList(i);
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
	m_edgeList.assign(m_NC, -1);
	for (int i = 0; i < NF; ++i)
	{
		if (m_faceList[i] != -1)
		{
			const std::vector<int>& faceEdge = topo.FaceEdgeList(i);
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

	feLog("\tRefinement info:\n");
	feLog("\t  Elements to refine: %d\n", m_splitElems);
	feLog("\t  Facets to refine  : %d\n", m_splitFaces);
	feLog("\t  Edges to refine   : %d\n", m_splitEdges);
}

int findNodeInMesh(FEMesh& mesh, const vec3d& r, double tol = 1e-12)
{
	int NN = mesh.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			vec3d ri = mesh.Node(i).m_r0;
			if ((ri - r).norm2() < tol) return i;
		}
	}
	return -1;
}

void FEHexRefine::UpdateNewNodes(FEModel& fem)
{
	FEMeshTopo& topo = *m_topo;
	FEMesh& mesh = fem.GetMesh();

	// we need to create a new node for each edge, face, and element that needs to be split
	int newNodes = m_splitEdges + m_splitFaces + m_splitElems;

	// for now, store the position of these new nodes in an array
	vector<vec3d> newPos(newNodes);

	// get the position of new nodes
	int n = 0;
	for (int i = 0; i < topo.Edges(); ++i)
	{
		if (m_edgeList[i] != -1)
		{
			const FEEdgeList::EDGE& edge = topo.Edge(i);
			vec3d r0 = mesh.Node(edge.node[0]).m_r0;
			vec3d r1 = mesh.Node(edge.node[1]).m_r0;
			newPos[n++] = (r0 + r1)*0.5;
		}
	}
	for (int i = 0; i < topo.Faces(); ++i)
	{
		if (m_faceList[i] != -1)
		{
			const FEFaceList::FACE& face = topo.Face(i);
			vec3d r0(0, 0, 0);
			int nn = face.ntype;
			for (int j = 0; j < nn; ++j) r0 += mesh.Node(face.node[j]).m_r0;
			r0 /= (double)nn;
			newPos[n++] = r0;
		}
	}
	for (int i = 0; i < topo.Elements(); ++i)
	{
		if (m_elemList[i] != -1)
		{
			FESolidElement& el = dynamic_cast<FESolidElement&>(*topo.Element(i));

			vec3d r0(0, 0, 0);
			int nn = el.Nodes();
			for (int j = 0; j < nn; ++j) r0 += mesh.Node(el.m_node[j]).m_r0;
			r0 /= (double)nn;
			newPos[n++] = r0;
		}
	}

	// some of these new nodes may coincide with an existing node
	//If we find one, we eliminate it
	n = 0;
	int nremoved = 0;
	for (int i = 0; i < topo.Edges(); ++i)
	{
		if (m_edgeList[i] != -1)
		{
			int nodeId = findNodeInMesh(mesh, newPos[n++]);
			if (nodeId >= 0)
			{
				if(mesh.Node(nodeId).HasFlags(FENode::HANGING))
					mesh.Node(nodeId).UnsetFlags(FENode::HANGING);
				m_edgeList[i] = -nodeId-2;
				nremoved++;
			}
		}
	}
	for (int i = 0; i < topo.Faces(); ++i)
	{
		if (m_faceList[i] != -1)
		{
			int nodeId = findNodeInMesh(mesh, newPos[n++]);
			if (nodeId >= 0)
			{
				if(mesh.Node(nodeId).HasFlags(FENode::HANGING))
					mesh.Node(nodeId).UnsetFlags(FENode::HANGING);
				m_faceList[i] = -nodeId-2;
				nremoved++;
			}
		}
	}
	for (int i = 0; i < topo.Elements(); ++i)
	{
		if (m_elemList[i] != -1)
		{
			int nodeId = findNodeInMesh(mesh, newPos[n++]);
			if (nodeId >= 0)
			{
				if(mesh.Node(nodeId).HasFlags(FENode::HANGING))
					mesh.Node(nodeId).UnsetFlags(FENode::HANGING);
				m_elemList[i] = -nodeId-2;
				nremoved++;
			}
		}
	}

	// we need to reindex nodes if some were removed
	if (nremoved > 0)
	{
		n = m_N0;
		for (int i = 0; i < topo.Edges(); ++i)
		{
			if (m_edgeList[i] >= 0) m_edgeList[i] = n++;
		}
		for (int i = 0; i < topo.Faces(); ++i)
		{
			if (m_faceList[i] >= 0) m_faceList[i] = n++;
		}
		for (int i = 0; i < topo.Elements(); ++i)
		{
			if (m_elemList[i] >= 0) m_elemList[i] = n++;
		}
		assert(n == (m_N0 + newNodes - nremoved));
		newNodes -= nremoved;
	}

	// now, generate new nodes
	mesh.AddNodes(newNodes);

	// assign dofs to new nodes
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	m_NN = mesh.Nodes();
	for (int i = m_N0; i<m_NN; ++i)
	{
		FENode& node = mesh.Node(i);
		node.SetDOFS(MAX_DOFS);
	}

	// update the position of these new nodes
	n = 0;
	for (int i = 0; i < topo.Edges(); ++i)
	{
		if (m_edgeList[i] >= 0)
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
		if (m_faceList[i] >= 0)
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
	for (int i=0; i < topo.Elements(); ++i)
	{
		if (m_elemList[i] >= 0)
		{
			FESolidElement& el = dynamic_cast<FESolidElement&>(*topo.Element(i));

			int nn = el.Nodes();
			FENode& node = mesh.Node(m_elemList[i]);

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
		if (m_edgeList[i] >= 0)
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
		if (m_faceList[i] >= 0)
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
	for (int i=0; i<topo.Elements(); ++i)
	{
		if (m_elemList[i] >= 0)
		{
			FESolidElement& el = dynamic_cast<FESolidElement&>(*topo.Element(i));
			int nn = el.Nodes();
			FENode& node = mesh.Node(m_elemList[i]);
			for (int j = 0; j < MAX_DOFS; ++j)
			{
				double v = 0.0;
				for (int k = 0; k < nn; ++k) v += mesh.Node(el.m_node[k]).get(j);
				v /= (double)nn;
				node.set(j, v);
			}
		}
	}

	// make the new node indices all positive
	for (int i = 0; i < topo.Edges(); ++i)
	{
		if (m_edgeList[i] < -1) m_edgeList[i] = -m_edgeList[i]-2;
	}
	for (int i = 0; i < topo.Faces(); ++i)
	{
		if (m_faceList[i] < -1) m_faceList[i] = -m_faceList[i]-2;
	}
	for (int i = 0; i < topo.Elements(); ++i)
	{
		if (m_elemList[i] < -1) m_elemList[i] = -m_elemList[i] - 2;
	}
}

//-----------------------------------------------------------------------------
// This function identifies the hanging nodes and assigns linear constraints to them.
void FEHexRefine::FindHangingNodes(FEModel& fem)
{
	FEMeshTopo& topo = *m_topo;
	FEMesh& mesh = fem.GetMesh();

	m_hangingNodes = 0;

	FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();

	const int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();

	// First, we removed any constraints on nodes that are no longer hanging
	int nremoved = 0;
	for (int i = 0; i < LCM.LinearConstraints();)
	{
		FELinearConstraint& lc = LCM.LinearConstraint(i);
		int nodeID = lc.master.node;
		if (mesh.Node(nodeID).HasFlags(FENode::HANGING) == false)
		{
			LCM.RemoveLinearConstraint(i);
			nremoved++;
		}
		else i++;
	}
	feLog("\tRemoved linear constraints : %d\n", nremoved);

	// we loop over non-split faces
	int nadded = 0;
	vector<int> tag(m_NC, 0);
	int NF = topo.Faces();
	for (int i = 0; i < NF; ++i)
	{
		if (m_faceList[i] == -1)
		{
			// This face is not split
			const std::vector<int>& fel = topo.FaceEdgeList(i);
			for (int j = 0; j < fel.size(); ++j)
			{
				if ((m_edgeList[fel[j]] >= 0) && (tag[fel[j]] == 0))
				{
					// Tag the node as hanging so we can identify it easier later
					int nodeId = m_edgeList[fel[j]];
					FENode& node = mesh.Node(nodeId);
					node.SetFlags(FENode::HANGING);

					// get the edge
					const FEEdgeList::EDGE& edge = topo.Edge(fel[j]);

					// setup a linear constraint for this node
					for (int k = 0; k < MAX_DOFS; ++k)
					{
						FELinearConstraint lc(&fem);
						lc.SetMasterDOF(k, nodeId);
						lc.AddSlaveDof(k, edge.node[0], 0.5);
						lc.AddSlaveDof(k, edge.node[1], 0.5);

						LCM.AddLinearConstraint(lc);
						nadded++;
					}

					// set a tag to avoid double-counting
					tag[fel[j]] = 1;

//					feLog("Added linear constraints to hanging node: %d (%d, %d)\n", nodeId, edge.node[0], edge.node[1]);

					m_hangingNodes++;
				}
			}
		}
	}
	// we loop over the non-split elements
	int NEL = topo.Elements();
	for (int i = 0; i<NEL; ++i)
	{
		if (m_elemList[i] == -1)
		{
			// This element is not split
			// If any of its faces are split, then the corresponding node
			// will be hanging
			const std::vector<int>& elface = topo.ElementFaceList(i);
			for (int j = 0; j < elface.size(); ++j)
			{
				if (m_faceList[elface[j]] >= 0)
				{
					// Tag the node as hanging so we can identify it easier later
					int nodeId = m_faceList[elface[j]];
					FENode& node = mesh.Node(nodeId);
					node.SetFlags(FENode::HANGING);

					// get the face
					const FEFaceList::FACE& face = topo.Face(elface[j]);

					// setup a linear constraint for this node
					for (int k = 0; k < MAX_DOFS; ++k)
					{
						FELinearConstraint lc(&fem);
						lc.SetMasterDOF(k, nodeId);
						lc.AddSlaveDof(k, face.node[0], 0.25);
						lc.AddSlaveDof(k, face.node[1], 0.25);
						lc.AddSlaveDof(k, face.node[2], 0.25);
						lc.AddSlaveDof(k, face.node[3], 0.25);

						LCM.AddLinearConstraint(lc);
						nadded++;
					}

//					feLog("Added linear constraints to hanging node: %d (%d, %d, %d, %d)\n", nodeId, face.node[0], face.node[1], face.node[2], face.node[3]);

					m_hangingNodes++;
				}
			}

			// This element is not split
			// If any of its edges are split, then the corresponding node
			// will be hanging
			const std::vector<int>& eledge = topo.ElementEdgeList(i);
			for (int j = 0; j < eledge.size(); ++j)
			{
				if ((m_edgeList[eledge[j]] >= 0) && (tag[eledge[j]] == 0))
				{
					// Tag the node as hanging so we can identify it easier later
					int nodeId = m_edgeList[eledge[j]];
					FENode& node = mesh.Node(nodeId);
					node.SetFlags(FENode::HANGING);

					// get the edge
					const FEEdgeList::EDGE& edge = topo.Edge(eledge[j]);

					// setup a linear constraint for this node
					for (int k = 0; k < MAX_DOFS; ++k)
					{
						FELinearConstraint lc(&fem);
						lc.SetMasterDOF(k, nodeId);
						lc.AddSlaveDof(k, edge.node[0], 0.5);
						lc.AddSlaveDof(k, edge.node[1], 0.5);

						LCM.AddLinearConstraint(lc);
						nadded++;
					}

					// set a tag to avoid double-counting
					tag[eledge[j]] = 1;

//					feLog("Added linear constraints to hanging node: %d (%d, %d)\n", nodeId, edge.node[0], edge.node[1]);

					m_hangingNodes++;
				}
			}

		}
	}

	feLog("\tHanging nodes ............ : %d\n", m_hangingNodes);
	feLog("\tAdded linear constraints . : %d\n", nadded);
}

void FEHexRefine::BuildNewDomains(FEModel& fem)
{
	// This lookup table defines how a hex will be split in eight smaller hexes
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
				el1.SetMatID(el0.GetMatID());
				el1.setStatus(el0.status());
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
					const std::vector<int>& ee = topo.ElementEdgeList(nelems); assert(ee.size() == 12);
					const std::vector<int>& ef = topo.ElementFaceList(nelems); assert(ef.size() == 6);

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
					ENL[ 8] = m_edgeList[ee[0]];
					ENL[ 9] = m_edgeList[ee[1]];
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

					// assign nodes to new elements
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

						el1.SetMatID(el0.GetMatID());
						el1.setStatus(el0.status());
					}
				}
				else
				{
					// if the element is not split, we just copy the nodes from
					// the old domain
					FEElement& el1 = oldDom.ElementRef(nel++);
					el1.SetMatID(el0.GetMatID());
					el1.setStatus(el0.status());
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
