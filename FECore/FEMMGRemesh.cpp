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
#include "FEMMGRemesh.h"
#include "FEModel.h"
#include "FEMeshTopo.h"
#include "FEMesh.h"
#include "FEDomain.h"
#include "FESolidDomain.h"
#include "FESurface.h"
#include "log.h"
#include "FEOctreeSearch.h"

#ifdef HAS_MMG
#include "mmg/mmg3d/libmmg3d.h"
bool build_mmg_mesh(MMG5_pMesh mmgMesg, MMG5_pSol mmgSol, FEMeshTopo& topo, FEMeshAdaptorCriterion* criterion, vector<double>& metric, double scale);
bool build_new_mesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, FEModel& fem, vector<double>& metric);
#endif

BEGIN_FECORE_CLASS(FEMMGRemesh, FERefineMesh)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_maxiter, "max_iters");
	ADD_PARAMETER(m_maxelem, "max_elems");
	ADD_PARAMETER(m_hmin, "min_element_size");
	ADD_PARAMETER(m_hausd, "hausdorff");
	ADD_PARAMETER(m_hgrad, "gradation");
	ADD_PROPERTY(m_criterion, "criterion", 0);
END_FECORE_CLASS();

FEMMGRemesh::FEMMGRemesh(FEModel* fem) : FERefineMesh(fem)
{
	m_maxiter = 1;
	m_maxelem = 0;

	m_scale = 0.5;
	m_hmin = 0.0;
	m_hausd = 0.01;
	m_hgrad = 1.3;

	m_criterion = nullptr;
}

bool FEMMGRemesh::Apply(int iteration)
{
	FEModel& fem = *GetFEModel();

	if (iteration >= m_maxiter)
	{
		feLogWarning("Max iterations reached.");
		return true;
	}

	// see if we should do anything
	FEMesh& mesh = fem.GetMesh();
	if ((m_maxelem > 0) && (mesh.Elements() >= m_maxelem))
	{
		feLog("\tElement limit reached.\n");
		return true;
	}

	// build the mesh-topo
	if (BuildMeshTopo() == false)
	{
		feLogError("Cannot apply tetgen refinement: Error building topo structure.");
		return true;
	}

	// do the mesh refinement
	if (Remesh() == false)
	{
		feLogError("Nothing to do.");
		return true;
	}

	// reactivate the model
	UpdateModel();

	// all done
	return false;
}

bool FEMMGRemesh::Remesh()
{
#ifdef HAS_MMG
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	if (mesh.IsType(ET_TET4) == false) return false;

	// initialize the MMG mesh
	MMG5_pMesh mmgMesh = NULL;
	MMG5_pSol  mmgSol = NULL;
	MMG3D_Init_mesh(MMG5_ARG_start,
		MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
		MMG5_ARG_end);

	// --- build the MMG mesh ---
	FEMeshTopo& topo = *m_topo;
	if (build_mmg_mesh(mmgMesh, mmgSol, topo, m_criterion, m_metric, m_scale) == false) return false;
	
	// set the control parameters
	MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmin, m_hmin);
	MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hausd, m_hausd);
	MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hgrad, m_hgrad);

	// run the mesher
	int ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);

	if (ier == MMG5_STRONGFAILURE) {
		feLogError("MMG was not able to remesh the mesh.");
		return false;
	}
	else if (ier == MMG5_LOWFAILURE)
	{
		feLogError("MMG return low failure error");
	}

	// build the new mesh
	bool bret = build_new_mesh(mmgMesh, mmgSol, fem, m_metric);

	// Clean up
	MMG3D_Free_all(MMG5_ARG_start,
		MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
		MMG5_ARG_end);

	return bret;

#else
	return false;
#endif
}

#ifdef HAS_MMG

bool build_mmg_mesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, FEMeshTopo& topo, FEMeshAdaptorCriterion* criterion, vector<double>& metric, double scale)
{
	FEMesh& mesh = *topo.GetMesh();
	int NN = mesh.Nodes();
	int NE = topo.Elements();
	int NF = topo.SurfaceFaces();

	// allocate mesh size
	if (MMG3D_Set_meshSize(mmgMesh, NN, NE, 0, NF, 0, 0) != 1)
	{
		assert(false);
		return nullptr;
	}

	// set the vertex coordinates
	for (int i = 0; i < NN; ++i)
	{
		FENode& vi = mesh.Node(i);
		vec3d r = vi.m_r0;
		MMG3D_Set_vertex(mmgMesh, r.x, r.y, r.z, 0, i + 1);
	}

	// set the tetrahedra
	int c = 1;
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int ne = dom.Elements();
		for (int j = 0; j < ne; ++j, ++c)
		{
			FEElement& e = dom.ElementRef(j);
			int* n = &e.m_node[0];
			MMG3D_Set_tetrahedron(mmgMesh, n[0] + 1, n[1] + 1, n[2] + 1, n[3] + 1, i, c);
		}
	}

	// set the facet markers
	vector<int> faceMarker(NF, 0);

	int faceMark = 1;
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		FESurface& surf = mesh.Surface(i);
		vector<int> faceIndexList = topo.SurfaceFaceIndexList(surf);
		for (int j = 0; j < faceIndexList.size(); ++j)
		{
			faceMarker[faceIndexList[j]] = faceMark;
		}
		faceMark++;
	}

	// for node sets we are going to create artificial surfaces
	vector<int> nodeTags;
	for (int i = 0; i < mesh.NodeSets(); ++i)
	{
		FENodeSet& nset = *mesh.NodeSet(i);
		if (nset.Size() != mesh.Nodes())
		{
			nodeTags.assign(mesh.Nodes(), 0);
			for (int j = 0; j < nset.Size(); ++j) nodeTags[nset[j]] = 2;

			// see if this is indeed a surface node set
			for (int j = 0; j < NF; ++j)
			{
				const FEFaceList::FACE& face = topo.SurfaceFace(j);
				const int* fn = face.node;
				if ((nodeTags[fn[0]] != 0) && (nodeTags[fn[1]] != 0) && (nodeTags[fn[2]] != 0))
				{
					nodeTags[fn[0]] = 1;
					nodeTags[fn[1]] = 1;
					nodeTags[fn[2]] = 1;
				}
			}
			int twos = 0;
			for (int j = 0; j < mesh.Nodes(); ++j) if (nodeTags[j] == 2) twos++;

			if (twos == 0)
			{
				for (int j = 0; j < NF; ++j)
				{
					const FEFaceList::FACE& face = topo.SurfaceFace(j);
					const int* fn = face.node;
					if ((nodeTags[fn[0]] == 1) && (nodeTags[fn[1]] == 1) && (nodeTags[fn[2]] == 1))
					{
						assert(faceMarker[j] == 0);
						faceMarker[j] = faceMark;
					}
				}
			}
			faceMark++;
		}
	}

	// create the faces
	for (int i = 0; i < NF; ++i)
	{
		const FEFaceList::FACE& f = topo.SurfaceFace(i);
		const int* n = &f.node[0];
		MMG3D_Set_triangle(mmgMesh, n[0] + 1, n[1] + 1, n[2] + 1, faceMarker[i], i + 1);
	}

	// Now, we build the "solution", i.e. the target element size.
	// If no elements are selected, we set a homogenous remeshing using the element size parameter.
	// set the "solution", i.e. desired element size
	if (MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, NN, MMG5_Scalar) != 1)
	{
		assert(false);
		return nullptr;
	}

	vector<double> edgeScale(NN, 1.0);
	if (metric.empty())
	{
		// build the edge length table
		int ET_TET[6][2] = { { 0,1 },{ 1,2 },{ 2,0 },{ 0,3 },{ 1,3 },{ 2,3 } };
		vector<pair<double, int> > edgeLength(NN, pair<double, int>(0.0, 0));
		for (int i = 0; i < NE; ++i)
		{
			FEElement& el = *topo.Element(i);
			for (int j = 0; j < 6; ++j)
			{
				int a = el.m_node[ET_TET[j][0]];
				int b = el.m_node[ET_TET[j][1]];

				vec3d ra = mesh.Node(a).m_r0;
				vec3d rb = mesh.Node(b).m_r0;

				double L = (ra - rb).norm2();

//				if (L > edgeLength[a].first) edgeLength[a].first = L;
//				if (L > edgeLength[b].first) edgeLength[b].first = L;

				edgeLength[a].first += L; edgeLength[a].second++;
				edgeLength[b].first += L; edgeLength[b].second++;
			}
		}

		metric.resize(NN, 0.0);
		for (int i = 0; i < NN; ++i)
		{
			if (edgeLength[i].second != 0)
			{
				edgeLength[i].first /= (double)edgeLength[i].second;
				edgeLength[i].first = sqrt(edgeLength[i].first);
			}
			metric[i] = edgeLength[i].first;
		}
	}

	// scale factors
	if (criterion)
	{
		vector< pair<int, double> > elemList = criterion->GetElementList();
		for (int i = 0; i < (int)elemList.size(); ++i)
		{
			FEElement& el = *topo.Element(elemList[i].first);
			for (int j = 0; j < el.Nodes(); ++j)
			{
				double s = scale;
				if (elemList[i].second != 0.0) s = elemList[i].second;
				edgeScale[el.m_node[j]] = s;
			}
		}
	}
	else
	{
		edgeScale.assign(NN, scale);
	}

	// set the new metric
	for (int k = 0; k < NN; k++) {
		MMG3D_Set_scalarSol(mmgSol, metric[k]*edgeScale[k], k + 1);
	}

	return true;
}

int findClosestNodes(FEMesh& mesh, const vec3d& x, int k, vector<int>& closestNodes)
{
	int N0 = mesh.Nodes();
	if (N0 < k) k = N0;

	vector<double> dist(k, 0.0);
	closestNodes.resize(k);
	int n = 0;
	for (int i = 0; i < N0; ++i)
	{
		vec3d ri = mesh.Node(i).m_r0 - x;
		double L2 = ri*ri;

		if (n == 0)
		{
			closestNodes[0] = i;
			dist[0] = L2;
			n++;
		}
		else if (L2 <= dist[n - 1])
		{
			int m;
			for (m = 0; m < n; ++m)
			{
				if (L2 <= dist[m])
				{
					break;
				}
			}

			if (n < k) n++;
			for (int l = n-1; l > m; l--)
			{
				closestNodes[l] = closestNodes[l - 1];
				dist[l] = dist[l - 1];
			}

			closestNodes[m] = i;
			dist[m] = L2;
		}
		else if (n < k)
		{
			closestNodes[n] = i;
			dist[n] = L2;
			n++;
		}
	}

	return n;
}

bool build_new_mesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, FEModel& fem, vector<double>& metric)
{
	FEMesh& mesh = fem.GetMesh();
	int N0 = mesh.Nodes();

	// get the new mesh sizes
	int nodes, elems, faces;
	MMG3D_Get_meshSize(mmgMesh, &nodes, &elems, NULL, &faces, NULL, NULL);

	// copy nodal positions
	vector<vec3d> nodePos0(nodes);
	for (int i = 0; i<nodes; ++i)
	{
		vec3d r;
		MMG3D_Get_vertex(mmgMesh, &r.x, &r.y, &r.z, NULL, NULL, NULL);
		nodePos0[i] = r;
	}

	// copy the metric
	metric.resize(nodes, 0.0);
	for (int i = 0; i < nodes; ++i)
	{
		MMG3D_Get_scalarSol(mmgSol, &metric[i]);
	}

	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	vector<vec3d> nodePos(nodes);
	vector<vector<double> > nodeVal(nodes, vector<double>(MAX_DOFS, 0.0));

	int transferMethod = 1;

	if (transferMethod == 0)
	{
		FEOctreeSearch octree(&mesh);
		octree.Init();

		// update solution
		for (int i = 0; i < nodes; ++i)
		{
			vec3d ri = nodePos0[i];

			double r[3] = { 0 };
			FESolidElement* el = (FESolidElement*)octree.FindElement(ri, r);
			if (el == nullptr)
			{
				assert(false);
				return false;
			}

			// get the nodal coordinates
			vec3d rt[FEElement::MAX_NODES];
			for (int j = 0; j < el->Nodes(); ++j) rt[j] = mesh.Node(el->m_node[j]).m_rt;

			nodePos[i] = el->evaluate(rt, r[0], r[1], r[2]);

			// update values
			for (int l = 0; l < MAX_DOFS; ++l)
			{
				double v[FEElement::MAX_NODES] = { 0 };
				for (int j = 0; j < el->Nodes(); ++j) v[j] = mesh.Node(el->m_node[j]).get(l);
				double vl = el->evaluate(v, r[0], r[1], r[2]);
				nodeVal[i][l] = vl;
			}
		}
	}
	else
	{
		// --- Do moving least-squares transfer ---

		// loop over all the new nodes
		for (int i = 0; i < nodes; ++i)
		{
			// get the spatial position
			vec3d x = nodePos0[i];

			// Find the closest M points
			vector<int> closestNodes;
			int M = findClosestNodes(mesh, x, 8, closestNodes);
			assert(M > 4);

			// the last node is the farthest and determines the radius
			vec3d& r = mesh.Node(closestNodes[M - 1]).m_r0;
			double L = sqrt((r - x)*(r - x));

			// add some offset to make sure none of the points will have a weight of zero.
			// Such points would otherwise be ignored, which reduces the net number of interpolation
			// points and make the MLS system ill-conditioned
			L += L*0.05;

			// evaluate weights and displacements
			vector<vec3d> Xi(M), Ui(M);
			vector<double> W(M);
			for (int m = 0; m < M; ++m)
			{
				FENode& node = mesh.Node(closestNodes[m]);
				vec3d uj = node.m_rt - node.m_r0;
				vec3d rj = x - node.m_r0;

				Xi[m] = rj;
				Ui[m] = uj;

				double D = sqrt(rj*rj);
				double wj = 1.0 - D / L;
				W[m] = wj;
			}

			// see if we have a (nearly) exact match with the closest point
			if (W[0] > 0.999)
			{
				FENode& node = mesh.Node(closestNodes[0]);

				// just use this point
				nodePos[i] = node.m_rt;
				for (int l = 0; l < MAX_DOFS; ++l)
				{
					nodeVal[i][l] = node.get(l);
				}
			}
			else
			{

				// setup least squares problems
				matrix A(4, 4);
				A.zero();
				vector<double> bx(4, 0.0), by(4, 0.0), bz(4, 0.0);
				for (int m = 0; m < M; ++m)
				{
					vec3d ri = Xi[m];
					vec3d ui = Ui[m];
					double P[4] = { 1.0, ri.x, ri.y, ri.z };

					for (int a = 0; a < 4; ++a)
					{
						for (int b = 0; b < 4; ++b)
						{
							A(a, b) += W[m] * P[a] * P[b];
						}

						bx[a] += W[m] * P[a] * ui.x;
						by[a] += W[m] * P[a] * ui.y;
						bz[a] += W[m] * P[a] * ui.z;
					}
				}

				// solve the linear system of equations
				vector<int> indx(4);
				A.lufactor(indx);
				A.lusolve(bx, indx);
				A.lusolve(by, indx);
				A.lusolve(bz, indx);

				// evaluate at x
				// this gives us the displacement
				vec3d ui(bx[0], by[0], bz[0]);

				// set the new position
				nodePos[i] = nodePos0[i] + ui;

				// update nodal values
				for (int l = 0; l < MAX_DOFS; ++l)
				{
					vector<double> b(4, 0.0);
					for (int m = 0; m < M; ++m)
					{
						vec3d ri = Xi[m];
						double P[4] = { 1.0, ri.x, ri.y, ri.z };

						double vm = mesh.Node(closestNodes[m]).get(l);

						for (int a = 0; a < 4; ++a)
						{
							b[a] += W[m] * P[a] * vm;
						}
					}

					A.lusolve(b, indx);

					nodeVal[i][l] = b[0];
				}
			}
		}
	}

	// reallocate nodes
	mesh.CreateNodes(nodes);

	// assign dofs to new nodes
	for (int i = 0; i < nodes; ++i)
	{
		FENode& node = mesh.Node(i);
		node.SetDOFS(MAX_DOFS);
		node.m_r0 = nodePos0[i];
		node.m_rt = nodePos[i];
		for (int j = 0; j < node.m_ID.size(); ++j) {
			node.set(j, nodeVal[i][j]);
		}
		node.UpdateValues();
	}

	// recreate domains
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(i));

		int nelems = 0;
		for (int j = 0; j < elems; ++j)
		{
			int n[4], gid, breq;
			MMG3D_Get_tetrahedron(mmgMesh, n, n + 1, n + 2, n + 3, &gid, &breq);
			if (gid == i) nelems++;
		}

		dom.Create(nelems, FE_TET4G1);
		int c = 0;
		for (int j = 0; j < elems; ++j)
		{
			int n[4], gid;
			MMG3D_Get_tetrahedron(mmgMesh, n, n + 1, n + 2, n + 3, &gid, NULL);

			if (gid == i)
			{
				FESolidElement& el = dom.Element(c++);
				el.m_node[0] = n[0] - 1;
				el.m_node[1] = n[1] - 1;
				el.m_node[2] = n[2] - 1;
				el.m_node[3] = n[3] - 1;
			}
		}
	}

	// re-init domains
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.CreateMaterialPointData();
		dom.Init();
		dom.Activate();
	}

	// recreate surfaces
	int faceMark = 1;
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		FESurface& surf = mesh.Surface(i);

		// count faces
		int nfaces = 0;
		for (int j = 0; j < faces; ++j)
		{
			int n[3], gid, breq;
			MMG3D_Get_triangle(mmgMesh, n, n + 1, n + 2, &gid, &breq);
			if (gid == faceMark) nfaces++;
		}
		assert(nfaces > 0);
		surf.Create(nfaces);
		int c = 0;
		for (int j = 0; j < faces; ++j)
		{
			int n[3], gid;
			MMG3D_Get_triangle(mmgMesh, n, n + 1, n + 2, &gid, NULL);
			if (gid == faceMark)
			{
				FESurfaceElement& face = surf.Element(c++);
				face.SetType(FE_TRI3G3);
				face.m_node[0] = n[0] - 1;
				face.m_node[1] = n[1] - 1;
				face.m_node[2] = n[2] - 1;
			}
		}
		assert(c == nfaces);
		surf.Init();

		faceMark++;
	}

	// update nodesets
	for (int i = 0; i < mesh.NodeSets(); ++i)
	{
		FENodeSet& nset = *mesh.NodeSet(i);
		if (nset.Size() != N0)
		{
			vector<int> nodeTags(mesh.Nodes(), 0);
			for (int j = 0; j < faces; ++j)
			{
				int n[3], gid;
				MMG3D_Get_triangle(mmgMesh, n, n + 1, n + 2, &gid, NULL);
				if (gid == faceMark)
				{
					nodeTags[n[0] - 1] = 1;
					nodeTags[n[1] - 1] = 1;
					nodeTags[n[2] - 1] = 1;
				}
			}

			std::vector<int> nodeList;
			for (int i = 0; i < nodeTags.size(); ++i) if (nodeTags[i] == 1) nodeList.push_back(i);

			if (nodeList.size() > 0)
			{
				nset.Clear();
				nset.Add(nodeList);
			}

			faceMark++;
		}
	}

	return true;
}
#endif
