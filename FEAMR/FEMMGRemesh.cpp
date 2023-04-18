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
#include "FEMMGRemesh.h"
#include <FECore/FEMeshTopo.h>
#include <FECore/FEModel.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FESurface.h>
#include <FECore/log.h>
#include <FECore/FEOctreeSearch.h>
#include <FECore/FENNQuery.h>
#include <FECore/FEMeshAdaptorCriterion.h>
#include <FECore/FEDomainMap.h>
#include "FELeastSquaresInterpolator.h"
#include "FEMeshShapeInterpolator.h"
#include "FEDomainShapeInterpolator.h"
#include <FECore/FECoreKernel.h>
#ifdef HAS_MMG
#include "mmg/mmg3d/libmmg3d.h"

class FEMMGRemesh::MMG
{
public:
	MMG(FEMMGRemesh* mmgRemesh) : m_mmgRemesh(mmgRemesh) {}
	bool build_mmg_mesh(MMG5_pMesh mmgMesg, MMG5_pSol mmgSol, FEMeshTopo& topo);
	bool build_new_mesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, FEModel& fem);

public:
	FEMMGRemesh*	m_mmgRemesh;
	std::vector<double>	m_metric;	// refinement metric
	std::vector<int>	m_nodeSetTag;	// surface tags for node sets
};

#endif

BEGIN_FECORE_CLASS(FEMMGRemesh, FERefineMesh)
	ADD_PARAMETER(m_hmin, "min_element_size");
	ADD_PARAMETER(m_hausd, "hausdorff");
	ADD_PARAMETER(m_hgrad, "gradation");
	ADD_PARAMETER(m_relativeSize, "relative_size");
	ADD_PARAMETER(m_meshCoarsen, "mesh_coarsen");
	ADD_PARAMETER(m_normalizeData, "normalize_data");
	ADD_PROPERTY(m_criterion, "criterion");
	ADD_PROPERTY(m_sfunc, "size_function", 0);
END_FECORE_CLASS();

FEMMGRemesh::FEMMGRemesh(FEModel* fem) : FERefineMesh(fem)
{
	m_maxelem = 0;
	m_relativeSize = true;
	m_meshCoarsen = false;
	m_normalizeData = false;

	m_hmin = 0.0;
	m_hausd = 0.01;
	m_hgrad = 1.3;

	m_criterion = nullptr;

	m_transferMethod = TRANSFER_MLQ;
	m_nnc = 8;

	m_sfunc = nullptr;

#ifdef HAS_MMG
	mmg = new FEMMGRemesh::MMG(this);
#endif
}

bool FEMMGRemesh::Init()
{
	FEMesh& mesh = GetMesh();
	if (mesh.IsType(ET_TET4) == false) return false;

	return FERefineMesh::Init();
}

bool FEMMGRemesh::RefineMesh()
{
#ifdef HAS_MMG
	FEMesh& mesh = GetMesh();

	// initialize the MMG mesh
	MMG5_pMesh mmgMesh = NULL;
	MMG5_pSol  mmgSol = NULL;
	MMG3D_Init_mesh(MMG5_ARG_start,
		MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
		MMG5_ARG_end);

	// --- build the MMG mesh ---
	FEMeshTopo& topo = *m_topo;
	if (mmg->build_mmg_mesh(mmgMesh, mmgSol, topo) == false) return false;
	
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
	bool bret = mmg->build_new_mesh(mmgMesh, mmgSol, *GetFEModel());

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

bool FEMMGRemesh::MMG::build_mmg_mesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, FEMeshTopo& topo)
{
	FEMesh& mesh = *topo.GetMesh();
	int NN = mesh.Nodes();
	int NE = topo.Elements();
	int NF = topo.SurfaceFaces();

	// allocate mesh size
	if (MMG3D_Set_meshSize(mmgMesh, NN, NE, 0, NF, 0, 0) != 1)
	{
		assert(false);
		return false;
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
	m_nodeSetTag.assign(mesh.NodeSets(), -1);
	for (int i = 0; i < mesh.NodeSets(); ++i)
	{
		FENodeSet& nset = *mesh.NodeSet(i);
		if (nset.Size() != mesh.Nodes())
		{
			vector<int> nodeTags(mesh.Nodes(), 0);
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
						if (faceMarker[j] == 0)
						{
							faceMarker[j] = faceMark;
							m_nodeSetTag[i] = faceMark;
						}
						else
						{
							if (m_nodeSetTag[i] == -1)
							{
								m_nodeSetTag[i] = faceMarker[j];
							}
							else if (faceMarker[j] != m_nodeSetTag[i])
							{
								return false;
							}
						}
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
		return false;
	}

	if (m_metric.empty())
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

		m_metric.resize(NN, 0.0);
		for (int i = 0; i < NN; ++i)
		{
			if (edgeLength[i].second != 0)
			{
				edgeLength[i].first /= (double)edgeLength[i].second;
				edgeLength[i].first = sqrt(edgeLength[i].first);
			}
			m_metric[i] = edgeLength[i].first;
		}
	}

	FEElementSet* elset = m_mmgRemesh->GetElementSet();
	if (elset)
	{
		// elements that are not in the element set will be flagged as required.
		FEElementIterator it(&mesh);
		int c = 1;
		for (; it.isValid(); ++it, ++c)
		{
			FEElement& el = *it;
			if (elset->Contains(el) == false)
			{
				MMG3D_Set_requiredTetrahedron(mmgMesh, c);
			}
		}
	}

	// scale factors
	vector<double> nodeScale(NN, 0.0);
	FEMeshAdaptorCriterion* criterion = m_mmgRemesh->GetCriterion();

	assert(criterion);
	if (criterion == nullptr) return false;
	FEMeshAdaptorSelection elemList = criterion->GetElementSelection(elset);

	// see if want to normalize the data
	if (m_mmgRemesh->m_normalizeData)
	{
		// Find data range
		double vmin, vmax;
		for (int i = 0; i < (int)elemList.size(); ++i)
		{
			double v = elemList[i].m_elemValue;
			if ((i == 0) || (v < vmin)) vmin = v;
			if ((i == 0) || (v > vmax)) vmax = v;
		}
		if (vmax == vmin) vmax++;

		// normalize data
		for (int i = 0; i < (int)elemList.size(); ++i)
		{
			double v = elemList[i].m_elemValue;
			elemList[i].m_elemValue = (v - vmin) / (vmax - vmin);
		}
	}

	// map to nodal data
	vector<int> tag(NN, 0);
	for (int i = 0; i < (int)elemList.size(); ++i)
	{
		FEElement& el = *mesh.FindElementFromID(elemList[i].m_elementId);
		for (int j = 0; j < el.Nodes(); ++j)
		{
			double s = elemList[i].m_elemValue;
			FEFunction1D* fs = m_mmgRemesh->m_sfunc;
			if (fs)
			{
				s = fs->value(s);
			}
			assert(s > 0.0);
			if (s <= 0.0) return false;
			nodeScale[el.m_node[j]] += s;
			tag[el.m_node[j]]++;
		}
	}
	for (int i = 0; i < NN; ++i)
	{
		if (tag[i] > 0) nodeScale[i] /= (double)tag[i];
		else nodeScale[i] = (m_mmgRemesh->m_relativeSize ? 1.0 : m_metric[i]);
	}

	// adjust for relative scale flag
	if (m_mmgRemesh->m_relativeSize)
	{
		for (int k = 0; k < NN; k++) {
			nodeScale[k] *= m_metric[k];
		}
	}

	// determine new size field
	bool meshCoarsen = m_mmgRemesh->m_meshCoarsen;
	for (int k = 0; k < NN; k++) 
	{
		double s = nodeScale[k];
		if ((meshCoarsen) || (s < m_metric[k])) m_metric[k] = s;
	}

	// set the new metric
	for (int k = 0; k < NN; k++) {
		MMG3D_Set_scalarSol(mmgSol, m_metric[k], k + 1);
	}

	return true;
}

bool FEMMGRemesh::MMG::build_new_mesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, FEModel& fem)
{
	FEMesh& mesh = fem.GetMesh();
	int N0 = mesh.Nodes();

	// get the new mesh sizes
	int nodes, elems, faces;
	MMG3D_Get_meshSize(mmgMesh, &nodes, &elems, NULL, &faces, NULL, NULL);

	// get old node positions
	vector<vec3d> oldNodePos(N0);
	for (int i = 0; i < N0; ++i) oldNodePos[i] = mesh.Node(i).m_r0;

	// copy nodal positions
	vector<vec3d> nodePos0(nodes);
	for (int i = 0; i<nodes; ++i)
	{
		vec3d r;
		MMG3D_Get_vertex(mmgMesh, &r.x, &r.y, &r.z, NULL, NULL, NULL);
		nodePos0[i] = r;
	}

	// copy the metric
	m_metric.resize(nodes, 0.0);
	for (int i = 0; i < nodes; ++i)
	{
		MMG3D_Get_scalarSol(mmgSol, &m_metric[i]);
	}

	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	vector<vec3d> nodePos(nodes);
	vector<vector<double> > nodeVal(nodes, vector<double>(MAX_DOFS, 0.0));

	// allocate data mapper
	FEMeshDataInterpolator* mapper = nullptr;
	switch (m_mmgRemesh->m_transferMethod)
	{
	case TRANSFER_SHAPE: mapper = new FEMeshShapeInterpolator(&mesh); break;
	case TRANSFER_MLQ:
	{
		FELeastSquaresInterpolator* MLQ = new FELeastSquaresInterpolator;
		MLQ->SetNearestNeighborCount(m_mmgRemesh->m_nnc);
		MLQ->SetDimension(m_mmgRemesh->m_nsdim);
		MLQ->SetSourcePoints(oldNodePos);
		mapper = MLQ;
	}
	break;
	default:
		assert(false);
		return false;
	}

	// map nodal positions and nodal data
	for (int i = 0; i < nodes; ++i)
	{
		vec3d ri = nodePos0[i];

		if (mapper->SetTargetPoint(ri) == false)
		{
			assert(false);
			delete mapper;
			return false;
		}

		// get the nodal coordinates
		nodePos[i] = mapper->MapVec3d([&mesh](int sourceNode) {
			return mesh.Node(sourceNode).m_rt;
		});

		// update values
		for (int l = 0; l < MAX_DOFS; ++l)
		{
			nodeVal[i][l] = mapper->Map([&mesh, l](int sourceNode) {
				return mesh.Node(sourceNode).get(l);
			});
		}
	}
	delete mapper;

	// reallocate nodes
	mesh.CreateNodes(nodes);

	// assign dofs to new nodes
	for (int i = 0; i < nodes; ++i)
	{
		FENode& node = mesh.Node(i);
		node.SetDOFS(MAX_DOFS);
		node.m_r0 = nodePos0[i];
		node.m_rt = nodePos[i];
		if (m_mmgRemesh->m_nsdim == 2) node.m_rt.z = node.m_r0.z;
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

		dom.Create(nelems, dom.GetElementSpec());
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

		// re-init domain
		dom.CreateMaterialPointData();
		dom.Reset();	// NOTE: we need to call this to actually call the Init function on the material points.
		dom.Init();
		dom.Activate();
	}
	mesh.RebuildLUT();

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

	// update nodesets
	for (int i = 0; i < mesh.NodeSets(); ++i)
	{
		int tag = m_nodeSetTag[i];
		FENodeSet& nset = *mesh.NodeSet(i);
		if (nset.Size() != N0)
		{
			vector<int> nodeTags(mesh.Nodes(), 0);
			for (int j = 0; j < faces; ++j)
			{
				int n[3], gid;
				MMG3D_Get_triangle(mmgMesh, n, n + 1, n + 2, &gid, NULL);
				if (gid == tag)
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
		else
		{
			// assume this node set is determined by the entire mesh
			FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(0));
//			if (nset.Size() == dom.Nodes())
			{
				nset.Clear();
				for (int i = 0; i < dom.Nodes(); ++i) nset.Add(dom.NodeIndex(i));
			}
		}
	}

	return true;
}

#endif
