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
#include "FEHexRefine2D.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEMeshTopo.h>
#include <FECore/FEPrescribedDOF.h>
#include <FECore/FEElementList.h>
#include <FECore/FELinearConstraint.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FEMeshAdaptorCriterion.h>
#include <FECore/FESurface.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEHexRefine2D, FERefineMesh)
	ADD_PARAMETER(m_elemRefine, "max_elem_refine");
	ADD_PROPERTY(m_criterion, "criterion");
	ADD_PARAMETER(m_maxValue, "max_value");
END_FECORE_CLASS();

FEHexRefine2D::FEHexRefine2D(FEModel* fem) : FERefineMesh(fem)
{
	m_maxValue = 0.0;
	m_elemRefine = 0;
	m_criterion = nullptr;
}

bool FEHexRefine2D::Init()
{
	FEMesh& mesh = GetMesh();

	if (mesh.IsType(ET_HEX8) == false)
	{
		feLogError("Cannot apply hex refinement: Mesh is not a HEX8 mesh.");
		return true;
	}

	return FERefineMesh::Init();
}

bool FEHexRefine2D::RefineMesh()
{
	FEMeshTopo& topo = *m_topo;

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = GetMesh();
	FEElementList allElems(mesh);

	const int NEL = mesh.Elements();
	const int NDOM = mesh.Domains();
	m_N0 = mesh.Nodes();
	m_NC = topo.Edges();
	int NF = topo.Faces();

	// Build the lists of items to split
	if (BuildSplitLists(fem) == false) return false;

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

	// update node sets
	for (int i = 0; i < mesh.NodeSets(); ++i)
	{
		FENodeSet& nset = *mesh.NodeSet(i);
		UpdateNodeSet(nset);
	}

	// update all surfaces
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		FESurface& surf = mesh.Surface(i);
		if (UpdateSurface(surf) == false) return false;
	}

	// update all element sets
	for (int i = 0; i < mesh.ElementSets(); ++i)
	{
		FEElementSet& set = mesh.ElementSet(i);
		if (UpdateElementSet(set) == false) return false;
	}

	return true;
}

bool FEHexRefine2D::BuildSplitLists(FEModel& fem)
{
	FEMeshTopo& topo = *m_topo;
	FEMesh& mesh = fem.GetMesh();

	// Get the elements that we need to refine
	int NEL = mesh.Elements();
	m_elemList.assign(NEL, -1);
	if (m_criterion)
	{
		FEMeshAdaptorSelection selection = m_criterion->GetElementSelection(GetElementSet());
		for (int i = 0; i < selection.size(); ++i)
		{
			if (selection[i].m_elemValue > m_maxValue)
			{
				int eid = selection[i].m_elementId;
				int lid = topo.GetElementIndexFromID(eid);
				m_elemList[lid] = 1;
			}
		}
	}
	else
	{
		// just do'em all
		m_elemList.assign(NEL, 1);
	}

	// We cannot split elements that have hanging nodes
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
	m_splitElems = 0;
	for (int i = 0; i < m_elemList.size(); ++i) {
		if (m_elemList[i] == 1) {
			m_splitElems++;
		}
	}
	if (m_splitElems == 0) return true;

	// make sure we don't exceed the max elements per refinement step
	if ((m_elemRefine > 0) && (m_splitElems > m_elemRefine))
	{
		m_splitElems = 0;
		for (int i = 0; i < m_elemList.size(); ++i) {
			if (m_elemList[i] == 0) {
				if (m_splitElems >= m_elemRefine)
				{
					m_elemList[i] = -1;
				}
				else
				{
					m_splitElems++;
				}
			}
		}
		assert(m_splitElems == m_elemRefine);
	}

	// figure out which faces to refine
	int NF = topo.Faces();
	m_faceList.assign(NF, -1);

	// Faces in the XY plane will be split in four. Other faces will be split in two. 
	// So, we need to figure out which faces lie in the XY plane and which do not. 
	// We will mark faces that are not in the XY plane by a -2
	for (int i = 0; i < NF; ++i)
	{
		const FEFaceList::FACE& face = topo.Face(i);

		// calculate face normal
		vec3d r0 = mesh.Node(face.node[0]).m_r0;
		vec3d r1 = mesh.Node(face.node[1]).m_r0;
		vec3d r2 = mesh.Node(face.node[2]).m_r0;
		vec3d n = (r1 - r0) ^ (r2 - r0); n.unit();
		if (fabs(n.z) < 0.999)
		{
			// this normal is not perpendicular to XY plane, so mark it
			m_faceList[i] = -2;
		}
	}

	for (int i = 0; i < m_elemList.size(); ++i)
	{
		if (m_elemList[i] != -1)
		{
			int splitFaces = 0;
			const std::vector<int>& elface = topo.ElementFaceList(i);
			for (int j = 0; j < elface.size(); ++j)
			{
				if (m_faceList[elface[j]] == -1)
				{
					// this face will be split in 4
					m_faceList[elface[j]] = 1;
					splitFaces++;
				}
				else if (m_faceList[elface[j]] == -2)
				{
					// this face will be split in 2
					m_faceList[elface[j]] = 2;
					splitFaces++;
				}
			}

			// There should always be at least two faces to split.
			if (splitFaces < 2)
			{
				feLog("Cannot refine element due to error.");
				return false;
			}
		}
	}

	// count how many faces to split
	int N1 = mesh.Nodes();
	m_splitFaces = 0;
	for (int i = 0; i < m_faceList.size(); ++i) 
	{
		if (m_faceList[i] == 1) {
			m_faceList[i] = N1++;
			m_splitFaces++;
		}
		else if (m_faceList[i] == 2)
		{
			m_splitFaces++;
			m_faceList[i] = -2;
		}
		else
		{
			m_faceList[i] = -1;
		}
	}

	// figure out which edges to refine
	m_edgeList.assign(m_NC, -1);
	for (int i = 0; i < NF; ++i)
	{
		if (m_faceList[i] >= 0)
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

	return true;
}

// in FEHexRefine.cpp
int findNodeInMesh(FEMesh& mesh, const vec3d& r, double tol = 1e-9);

void FEHexRefine2D::UpdateNewNodes(FEModel& fem)
{
	FEMeshTopo& topo = *m_topo;
	FEMesh& mesh = fem.GetMesh();

	// we need to create a new node for each edge and face that needs to be split
	int newNodes = m_splitEdges + m_splitFaces;

	// for now, store the position of these new nodes in an array
	vector<vec3d> newPos(newNodes);

	// get the position of new nodes
	int n = 0;
	for (int i = 0; i < topo.Edges(); ++i)
	{
		if (m_edgeList[i] >= 0)
		{
			const FEEdgeList::EDGE& edge = topo.Edge(i);
			vec3d r0 = mesh.Node(edge.node[0]).m_r0;
			vec3d r1 = mesh.Node(edge.node[1]).m_r0;
			newPos[n++] = (r0 + r1)*0.5;
		}
	}
	for (int i = 0; i < topo.Faces(); ++i)
	{
		if (m_faceList[i] >= 0)
		{
			const FEFaceList::FACE& face = topo.Face(i);
			vec3d r0(0, 0, 0);
			int nn = face.ntype;
			for (int j = 0; j < nn; ++j) r0 += mesh.Node(face.node[j]).m_r0;
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
		if (m_edgeList[i] >= 0)
		{
			int nodeId = findNodeInMesh(mesh, newPos[n++]);
			if (nodeId >= 0)
			{
				assert(nodeId < m_N0);
				if(mesh.Node(nodeId).HasFlags(FENode::HANGING))
					mesh.Node(nodeId).UnsetFlags(FENode::HANGING);
				m_edgeList[i] = nodeId;
				nremoved++;
			}
		}
	}
	for (int i = 0; i < topo.Faces(); ++i)
	{
		if (m_faceList[i] >= 0)
		{
			int nodeId = findNodeInMesh(mesh, newPos[n++]);
			if (nodeId >= 0)
			{
				assert(nodeId < m_N0);
				if(mesh.Node(nodeId).HasFlags(FENode::HANGING))
					mesh.Node(nodeId).UnsetFlags(FENode::HANGING);
				m_faceList[i] = nodeId;
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
			if (m_edgeList[i] >= m_N0) m_edgeList[i] = n++;
		}
		for (int i = 0; i < topo.Faces(); ++i)
		{
			if (m_faceList[i] >= m_N0) m_faceList[i] = n++;
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
		if (m_edgeList[i] >= m_N0)
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
		if (m_faceList[i] >= m_N0)
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

	// re-evaluate solution at nodes
	for (int i = 0; i < topo.Edges(); ++i)
	{
		if (m_edgeList[i] >= m_N0)
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
		if (m_faceList[i] >= m_N0)
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
}

//-----------------------------------------------------------------------------
// This function identifies the hanging nodes and assigns linear constraints to them.
void FEHexRefine2D::FindHangingNodes(FEModel& fem)
{
	FEMeshTopo& topo = *m_topo;
	FEMesh& mesh = fem.GetMesh();

	m_hangingNodes = 0;

	FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();

	// First, we remove any constraints on nodes that are no longer hanging
	int nremoved = 0;
	for (int i = 0; i < LCM.LinearConstraints();)
	{
		FELinearConstraint& lc = LCM.LinearConstraint(i);
		int nodeID = lc.GetParentNode();
		if (mesh.Node(nodeID).HasFlags(FENode::HANGING) == false)
		{
			LCM.RemoveLinearConstraint(i);
			nremoved++;
		}
		else i++;
	}
	feLog("\tRemoved linear constraints : %d\n", nremoved);

	// Total nr of degrees of freedom
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();

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
						FELinearConstraint* lc = new FELinearConstraint(&fem);
						lc->SetParentDof(k, nodeId);
						lc->AddChildDof(k, edge.node[0], 0.5);
						lc->AddChildDof(k, edge.node[1], 0.5);

						LCM.AddLinearConstraint(lc);
						nadded++;
					}

					// set a tag to avoid double-counting
					tag[fel[j]] = 1;

					m_hangingNodes++;
				}
			}
		}
	}

	feLog("\tHanging nodes ............ : %d\n", m_hangingNodes);
	feLog("\tAdded linear constraints . : %d\n", nadded);
}

void FEHexRefine2D::BuildNewDomains(FEModel& fem)
{
	// This lookup table defines how a hex will be split in four smaller hexes
	const int LUT[4][8] = {
		{  0,  8, 16, 11,  4, 12, 17, 15 },
		{  8,  1,  9, 16, 12,  5, 13, 17 },
		{ 11, 16, 10,  3, 15, 17, 14,  7 },
		{ 16,  9,  2, 10, 17, 13,  6, 14 }
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
			FEDomain* newDom = fecore_new<FESolidDomain>(oldDom.GetTypeStr(), &fem);
			newDom->Create(NE0, FEElementLibrary::GetElementSpecFromType(FE_HEX8G8));
			for (int j = 0; j < NE0; ++j)
			{
				FEElement& el0 = oldDom.ElementRef(j);
				FEElement& el1 = newDom->ElementRef(j);
				for (int k = 0; k < el0.Nodes(); ++k) el1.m_node[k] = el0.m_node[k];
				el1.m_val = el0.m_val;
			}

			// reallocate the old domain
			oldDom.Create(4 * newElems + (NE0 - newElems), FEElementLibrary::GetElementSpecFromType(FE_HEX8G8));

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
					ENL[16] = m_faceList[ef[4]]; assert(ENL[16] >= 0);
					ENL[17] = m_faceList[ef[5]]; assert(ENL[17] >= 0);

					// assign nodes to new elements
					for (int k = 0; k < 4; ++k)
					{
						FEElement& el1 = oldDom.ElementRef(nel++);
						el1.m_val = el0.m_val;

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
					// if the element is not split, we just copy the nodes from
					// the old domain
					FEElement& el1 = oldDom.ElementRef(nel++);
					for (int k = 0; k < el0.Nodes(); ++k) el1.m_node[k] = el0.m_node[k];
					el1.m_val = el0.m_val;
				}
			}

			// we don't need this anymore
			delete newDom;
		}
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
}

void FEHexRefine2D::UpdateNodeSet(FENodeSet& nset)
{
	FEMeshTopo& topo = *m_topo;
	vector<int> tag(m_NN, 0);

	for (int j = 0; j < nset.Size(); ++j) tag[nset[j]] = 1;

	for (int j = 0; j < topo.Edges(); ++j)
	{
		if (m_edgeList[j] >= 0)
		{
			const FEEdgeList::EDGE& edge = topo.Edge(j);

			if ((tag[edge.node[0]] == 1) && (tag[edge.node[1]] == 1))
			{
				nset.Add(m_edgeList[j]);
			}
		}
	}

	for (int j = 0; j < topo.Faces(); ++j)
	{
		if (m_faceList[j] >= 0)
		{
			const FEFaceList::FACE& face = topo.Face(j);

			assert(face.ntype == 4);
			if ((tag[face.node[0]] == 1) &&
				(tag[face.node[1]] == 1) &&
				(tag[face.node[2]] == 1) &&
				(tag[face.node[3]] == 1))
			{
				nset.Add(m_faceList[j]);
			}
		}
	}
}

bool FEHexRefine2D::UpdateSurface(FESurface& surf)
{
	// look-up table for splitting quads in 4
	const int LUT[4][4] = {
		{ 0, 4, 8, 7 },
		{ 4, 1, 5, 8 },
		{ 7, 8, 6, 3 },
		{ 8, 5, 2, 6 }
	};

	FEMeshTopo& topo = *m_topo;
	int NF0 = surf.Elements();

	// figure out which facets to split
	vector<int> faceList = topo.FaceIndexList(surf);
	assert((int)faceList.size() == NF0);

	// count how many faces to split
	int split4 = 0, split2 = 0;
	for (int i = 0; i < faceList.size(); ++i)
	{
		int iface = faceList[i];
		if (m_faceList[iface] >= 0) split4++;
		if (m_faceList[iface] == -2) split2++;
	}
	if (split4 + split2 == 0) return surf.Init();

	// create a copy of the domain
	FESurface oldSurf(GetFEModel());
	oldSurf.Create(NF0);
	for (int i = 0; i < NF0; ++i)
	{
		FESurfaceElement& el0 = surf.Element(i);
		FESurfaceElement& el1 = oldSurf.Element(i);
		el1.SetType(el0.Type());
		int nf = el0.Nodes();
		for (int j = 0; j < nf; ++j) el1.m_node[j] = el0.m_node[j];
	}

	// reallocate the domain (Assumes Quad faces!)
	int NF1 = NF0 - split4 + 4 * (split4) - split2 + 2*split2;
	surf.Create(NF1);

	// reinitialize the surface
	int n = 0;
	for (int i = 0; i < NF0; ++i)
	{
		FESurfaceElement& el0 = oldSurf.Element(i);

		int iface = faceList[i];
		if (m_faceList[iface] >= 0)
		{
			const FEFaceList::FACE& face = topo.Face(iface);
			const vector<int>& edge = topo.FaceEdgeList(iface);

			int NL[9];
			NL[0] = face.node[0];
			NL[1] = face.node[1];
			NL[2] = face.node[2];
			NL[3] = face.node[3];
			NL[4] = m_edgeList[edge[0]];
			NL[5] = m_edgeList[edge[1]];
			NL[6] = m_edgeList[edge[2]];
			NL[7] = m_edgeList[edge[3]];
			NL[8] = m_faceList[iface];

			for (int j = 0; j < 4; ++j)
			{
				FESurfaceElement& el1 = surf.Element(n++);
				el1.SetType(FE_QUAD4G4);
				el1.m_node[0] = NL[LUT[j][0]];
				el1.m_node[1] = NL[LUT[j][1]];
				el1.m_node[2] = NL[LUT[j][2]];
				el1.m_node[3] = NL[LUT[j][3]];
			}
		}
		else if (m_faceList[iface] == -2)
		{
			const FEFaceList::FACE& face = topo.Face(iface);
			const vector<int>& edge = topo.FaceEdgeList(iface);

			int NL[2][4];

			// there should be two edges that are split, and two that are not
			int eid[4] = { m_edgeList[edge[0]], m_edgeList[edge[1]], m_edgeList[edge[2]], m_edgeList[edge[3]] };

			if ((eid[0] >= 0) && (eid[2] >= 0))
			{
				assert((eid[1] == -1) && (eid[3] == -1));
				NL[0][0] = face.node[0]; NL[1][0] = face.node[1];
				NL[0][1] =       eid[0]; NL[1][1] = face.node[2];
				NL[0][2] =       eid[2]; NL[1][2] =       eid[2];
				NL[0][3] = face.node[3]; NL[1][3] =       eid[0];
			}
			else if ((eid[1] >= 0) && (eid[3] >= 0))
			{
				assert((eid[0] == -1) && (eid[2] == -1));
				NL[0][0] = face.node[0]; NL[1][0] = face.node[2];
				NL[0][1] = face.node[1]; NL[1][1] = face.node[3];
				NL[0][2] =       eid[1]; NL[1][2] =       eid[3];
				NL[0][3] =       eid[3]; NL[1][3] =       eid[1];
			}
			else { assert(false); }

			for (int j = 0; j < 2; ++j)
			{
				FESurfaceElement& el1 = surf.Element(n++);
				el1.SetType(FE_QUAD4G4);
				el1.m_node[0] = NL[j][0];
				el1.m_node[1] = NL[j][1];
				el1.m_node[2] = NL[j][2];
				el1.m_node[3] = NL[j][3];
			}
		}
		else
		{
			FESurfaceElement& el1 = surf.Element(n++);
			el1.SetType(FE_QUAD4G4);
			el1.m_node[0] = el0.m_node[0];
			el1.m_node[1] = el0.m_node[1];
			el1.m_node[2] = el0.m_node[2];
			el1.m_node[3] = el0.m_node[3];
		}
	}
	surf.CreateMaterialPointData();

	return surf.Init();
}

bool FEHexRefine2D::UpdateElementSet(FEElementSet& eset)
{
	// get the domain list
	// NOTE: Don't get the reference, since then the same reference
	// is passed to Create below, which causes problems.
	FEDomainList domList = eset.GetDomainList();
	if (domList.IsEmpty()) { throw std::runtime_error("Error in FEHexRefine2D!"); }

	// recreate the element set from the domain list
	eset.Create(domList);

	return true;
}
