/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include <FECore/FEModel.h>
#include <FECore/FEMeshTopo.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FESurface.h>
#include <FECore/log.h>
#include <FECore/FEOctreeSearch.h>
#include <FECore/FENNQuery.h>
#include <FECore/FEMeshAdaptorCriterion.h>
#include <FECore/DumpMemStream.h>
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
	MMG(FEMMGRemesh* mmgRemesh) : m_mmgRemesh(mmgRemesh) 
	{
		m_meshCopy = nullptr;
	}
	~MMG()
	{
		clear_maps();
	}

	bool build_mmg_mesh(MMG5_pMesh mmgMesg, MMG5_pSol mmgSol, FEMeshTopo& topo, double scale);
	bool build_new_mesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, FEModel& fem);

private:
	bool build_map_data(FEModel& fem);
	void map_data(FEModel& fem);
	FEDomainMap* createElemDataMap(FEModel& fem, FEDomain& dom, vector<vec3d>& nodePos, FEDomainMap* map, FEMeshDataInterpolator* dataMapper);
	void NodeToElemData(FEModel& fem, FEDomain& dom, FEDomainMap* nodeMap, FEDomainMap* elemMap, FEMeshDataInterpolator* dataMapper);

	void clear_maps()
	{
		for (size_t i = 0; i < m_nodeMapList.size(); ++i)
		{
			std::vector<FEDomainMap*>& map_i = m_nodeMapList[i];
			for (size_t j = 0; j < map_i.size(); ++j) delete map_i[j];
		}
		m_nodeMapList.clear();

		for (int i = 0; i < m_meshDataList.size(); ++i) delete m_meshDataList[i];
		m_meshDataList.clear();
	}

public:
	FEMesh*		m_meshCopy;				// (shallow) copy of old mesh
	FEMMGRemesh*	m_mmgRemesh;
	std::vector<double>	m_metric;	// refinement metric
	std::vector<int>	m_nodeSetTag;	// surface tags for node sets

	std::vector< std::vector<vec3d> >		m_oldNodePos;	// old domain's nodal coordinates
	std::vector< std::vector<FEDomainMap*> >	m_nodeMapList;	// list of nodal data for each domain
	std::vector< FEDomainMap* >	m_meshDataList;					// list of nodal data for mesh data
};

#endif

BEGIN_FECORE_CLASS(FEMMGRemesh, FERefineMesh)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_maxiter, "max_iters");
	ADD_PARAMETER(m_maxelem, "max_elems");
	ADD_PARAMETER(m_hmin, "min_element_size");
	ADD_PARAMETER(m_hausd, "hausdorff");
	ADD_PARAMETER(m_hgrad, "gradation");
	ADD_PARAMETER(m_bmap_data, "map_data");
	ADD_PARAMETER(m_nnc, "nnc");
	ADD_PARAMETER(m_transferMethod, "transfer_method");
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

	m_transferMethod = TRANSFER_MLQ;

	m_bmap_data = false;
	m_nnc = 8;

#ifdef HAS_MMG
	mmg = new FEMMGRemesh::MMG(this);
#endif
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

	// print some mesh statistics
	int NN = mesh.Nodes();
	int NE = mesh.Elements();
	feLog(" Mesh Statistics:\n");
	feLog(" \tNumber of nodes    : %d\n", NN);
	feLog(" \tNumber of elements : %d\n", NE);
	feLog("\n");

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
	if (mmg->build_mmg_mesh(mmgMesh, mmgSol, topo, m_scale) == false) return false;
	
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
	bool bret = mmg->build_new_mesh(mmgMesh, mmgSol, fem);

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

bool FEMMGRemesh::MMG::build_mmg_mesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, FEMeshTopo& topo, double scale)
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

	vector<double> edgeScale(NN, 1.0);
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
	FEMeshAdaptorCriterion* criterion = m_mmgRemesh->GetCriterion();
	if (criterion)
	{
		FEMeshAdaptorSelection elemList = criterion->GetElementSelection(elset);
		for (int i = 0; i < (int)elemList.size(); ++i)
		{
			FEElement& el = *topo.Element(elemList[i].m_elementIndex);
			for (int j = 0; j < el.Nodes(); ++j)
			{
				double s = scale;
				if (elemList[i].m_scaleFactor != 0.0) s = elemList[i].m_scaleFactor;
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
		MMG3D_Set_scalarSol(mmgSol, m_metric[k]*edgeScale[k], k + 1);
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

	// before we recreate the mesh, we need to extract the data that
	// needs to be mapped between meshes.
	if (m_mmgRemesh->m_bmap_data)
	{
		if (build_map_data(fem) == false)
		{
			return false;
		}
	}

	// --- create a copy of the old mesh --- 
	m_meshCopy = new FEMesh(nullptr);
	m_meshCopy->CreateNodes(N0);
	for (int i = 0; i < N0; ++i)
	{
		m_meshCopy->Node(i) = mesh.Node(i);
	}

	// now allocate domains
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		const char* sz = dom.GetTypeStr();

		// create a new domain
		FEDomain* pd = fecore_new<FEDomain>(sz, nullptr);
		assert(pd);
		pd->SetMesh(m_meshCopy);

		// copy domain data
		pd->CopyFrom(&dom);

		// add it to the mesh
		m_meshCopy->AddDomain(pd);
	}
	m_meshCopy->RebuildLUT();

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
		MLQ->SetCheckForMatch(true);
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

		dom.Create(nelems, FEElementLibrary::GetElementSpecFromType(FE_TET4G1));
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
		dom.Init();
		dom.Activate();
	}
	mesh.RebuildLUT();

	// map data onto new mesh
	if (m_mmgRemesh->m_bmap_data)
	{
		map_data(fem);
	}

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

FEDomainMap* createNodeDataMap(FEDomain& dom, FEDomainMap* map)
{
	FEDataType dataType = map->DataType();
	int dataSize = 0;
	switch (dataType)
	{
	case FEDataType::FE_DOUBLE: dataSize = 1; break;
	case FEDataType::FE_VEC3D : dataSize = 3; break;
	case FEDataType::FE_MAT3D : dataSize = 9; break;
	case FEDataType::FE_MAT3DS: dataSize = 6; break;
	default:
		assert(false);
		return nullptr;
	}

	// temp storage 
	double si[FEElement::MAX_INTPOINTS*9];
	double sn[FEElement::MAX_NODES*9];

	// allocate node data
	int NN = dom.Nodes();
	vector<double> nodeData(NN*dataSize);

	// build tag list
	vector<int> tag(NN, 0);
	int NE = dom.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int ne = e.Nodes();
		for (int k = 0; k < ne; ++k)
		{
			tag[e.m_lnode[k]]++;
		}
	}

	// get the data format
	int dataFormat = map->StorageFormat();
	if ((dataFormat != FMT_MATPOINTS) && (dataFormat != FMT_MULT)) return nullptr;

	// loop over all elements
	for (int i = 0; i < NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int ne = e.Nodes();

		int ni = (dataFormat == FMT_MATPOINTS ? e.GaussPoints() : ne);

		for (int j = 0; j < dataSize; ++j)
		{
			// get the integration point values
			for (int k = 0; k < ni; ++k)
			{
				switch (dataType)
				{
				case FEDataType::FE_DOUBLE:
					si[k] = map->value<double>(i, k);
					break;
				case FEDataType::FE_VEC3D:
					{
						vec3d v = map->value<vec3d>(i, k);
						if (j == 0) si[k] = v.x;
						if (j == 1) si[k] = v.y;
						if (j == 2) si[k] = v.z;
					}
					break;
				case FEDataType::FE_MAT3D:
					{
						mat3d v = map->value<mat3d>(i, k);
						int LUT[9][2] = { {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {1,2}, {2,0}, {2,1}, {2,2} };
						si[k] = v(LUT[j][0], LUT[j][1]);
					}
					break;
				case FEDataType::FE_MAT3DS:
					{
						mat3ds v = map->value<mat3ds>(i, k);
						int LUT[6][2] = { {0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2} };
						si[k] = v(LUT[j][0], LUT[j][1]);
					}
					break;
				}
			}

			// project to nodes
			if (dataFormat == FMT_MATPOINTS)
			{
				e.project_to_nodes(si, sn);
			}
			else
			{
				for (int k = 0; k < ne; ++k) sn[k] = si[k];
			}

			for (int k = 0; k < ne; ++k)
			{
				nodeData[e.m_lnode[k]*dataSize + j] += sn[k];
			}
		}
	}

	// normalize data
	for (int i = 0; i < NN; ++i)
	{
		if (tag[i] > 0)
		{
			for (int j=0; j<dataSize; ++j)
				nodeData[i*dataSize + j] /= (double)tag[i];
		}
	}

	// create new data map
	FEDomainMap* nodeMap = new FEDomainMap(dataType, Storage_Fmt::FMT_NODE);
	nodeMap->Create(const_cast<FEElementSet*>(map->GetElementSet()));
	for (int i = 0; i < NN; ++i)
	{
		switch (dataType)
		{
		case FEDataType::FE_DOUBLE: nodeMap->setValue(i, nodeData[i]); break;
		case FEDataType::FE_VEC3D:
		{
			vec3d v;
			v.x = nodeData[i*dataSize  ];
			v.y = nodeData[i*dataSize+1];
			v.z = nodeData[i*dataSize+2];
			nodeMap->setValue(i, v);
		}
		break;
		case FEDataType::FE_MAT3D:
		{
			mat3d v;
			v(0, 0) = nodeData[i*dataSize    ]; v(0, 1) = nodeData[i*dataSize + 1]; v(0, 2) = nodeData[i*dataSize + 2];
			v(1, 0) = nodeData[i*dataSize + 3]; v(1, 1) = nodeData[i*dataSize + 4]; v(1, 2) = nodeData[i*dataSize + 5];
			v(2, 0) = nodeData[i*dataSize + 6]; v(2, 1) = nodeData[i*dataSize + 7]; v(2, 2) = nodeData[i*dataSize + 8];
			nodeMap->setValue(i, v);
		}
		break;
		case FEDataType::FE_MAT3DS:
		{
			mat3ds v;
			v(0, 0) = nodeData[i*dataSize    ]; 
			v(0, 1) = nodeData[i*dataSize + 1]; 
			v(0, 2) = nodeData[i*dataSize + 2];
			v(1, 1) = nodeData[i*dataSize + 3];
			v(1, 2) = nodeData[i*dataSize + 4];
			v(2, 2) = nodeData[i*dataSize + 5];
			nodeMap->setValue(i, v);
		}
		break;
		}
	}

	return nodeMap;
}

FEDomainMap* FEMMGRemesh::MMG::createElemDataMap(FEModel& fem, FEDomain& dom, vector<vec3d>& nodePos, FEDomainMap* map, FEMeshDataInterpolator* dataMapper)
{
	assert(map->StorageFormat() == Storage_Fmt::FMT_NODE);

	FEDataType dataType = map->DataType();
	int dataSize = 0;
	switch (dataType)
	{
	case FEDataType::FE_DOUBLE: dataSize = 1; break;
	case FEDataType::FE_VEC3D: dataSize = 3; break;
	case FEDataType::FE_MAT3D: dataSize = 9; break;
	case FEDataType::FE_MAT3DS: dataSize = 6; break;
	default:
		assert(false);
		return nullptr;
	}

	// count nr of integration points
	int NMP = 0;
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		NMP += el.GaussPoints();
	}

	int N0 = nodePos.size();

	// create new domain map
	FEDomainMap* elemData = new FEDomainMap(map->DataType(), Storage_Fmt::FMT_MATPOINTS);
	FEElementSet* eset = new FEElementSet(&fem);
	eset->Create(&dom);
	elemData->Create(eset);

	vector<double> srcData(N0);
	vector<double> trgData(NMP);

	vector< vector<double> > mappedData(NMP, vector<double>(9, 0.0));

	// loop over all the new nodes
	for (int l = 0; l < dataSize; ++l)
	{
		for (int i = 0; i < N0; ++i)
		{
			double vm = 0.0;
			switch (dataType)
			{
			case FEDataType::FE_DOUBLE: vm = map->value<double>(0, i); break;
			case FEDataType::FE_VEC3D:
			{
				vec3d v = map->value<vec3d>(0, i);
				if (l == 0) vm = v.x;
				if (l == 1) vm = v.y;
				if (l == 2) vm = v.z;
			}
			break;
			case FEDataType::FE_MAT3D:
			{
				int LUT[9][2] = { {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {1,2}, {2,0}, {2,1}, {2,2} };
				mat3d v = map->value<mat3d>(0, i);
				vm = v(LUT[l][0], LUT[l][1]);
			}
			break;
			case FEDataType::FE_MAT3DS:
			{
				int LUT[6][2] = { {0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2} };
				mat3ds v = map->value<mat3ds>(0, i);
				vm = v(LUT[l][0], LUT[l][1]);
			}
			break;
			default:
				assert(false);
			}
			srcData[i] = vm;
		}

		dataMapper->Map(trgData, [&srcData](int sourcePoint) {
			return srcData[sourcePoint];
		});

		for (int i = 0; i < NMP; ++i)
		{
			mappedData[i][l] = trgData[i];
		}
	}

	// write mapped data to domain map
	int n = 0;
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		for (int j = 0; j < nint; ++j)
		{
			vector<double>& vj = mappedData[n++];

			for (int l=0; l<dataSize;++l)
			{
				switch (dataType)
				{
				case FEDataType::FE_DOUBLE: elemData->setValue(i, j, vj[0]); break;
				case FEDataType::FE_VEC3D:
				{
					vec3d v;
					v.x = vj[0];
					v.y = vj[1];
					v.z = vj[2];
					elemData->setValue(i, j, v);
				}
				break;
				case FEDataType::FE_MAT3D:
				{
					mat3d v;
					v(0, 0) = vj[0]; v(0, 1) = vj[1]; v(0, 2) = vj[2];
					v(1, 0) = vj[3]; v(1, 1) = vj[4]; v(1, 2) = vj[5];
					v(2, 0) = vj[6]; v(2, 1) = vj[7]; v(2, 2) = vj[8];
					elemData->setValue(i, j, v);
				}
				break;
				case FEDataType::FE_MAT3DS:
				{
					mat3ds v;
					v(0, 0) = vj[0];
					v(0, 1) = vj[1];
					v(0, 2) = vj[2];
					v(1, 1) = vj[3];
					v(1, 2) = vj[4];
					v(2, 2) = vj[5];
					elemData->setValue(i, j, v);
				}
				break;
				default:
					assert(false);
				}
			}
		}
	}

	return elemData;
}

void FEMMGRemesh::MMG::NodeToElemData(FEModel& fem, FEDomain& dom, FEDomainMap* nodeMap, FEDomainMap* elemMap, FEMeshDataInterpolator* dataMapper)
{
	assert(nodeMap->StorageFormat() == Storage_Fmt::FMT_NODE);

	FEDataType dataType = nodeMap->DataType();
	int dataSize = 0;
	switch (dataType)
	{
	case FEDataType::FE_DOUBLE: dataSize = 1; break;
	case FEDataType::FE_VEC3D: dataSize = 3; break;
	case FEDataType::FE_MAT3D: dataSize = 9; break;
	case FEDataType::FE_MAT3DS: dataSize = 6; break;
	default:
		assert(false);
		throw std::runtime_error("Error in FEMMGRemesh::MMG::NodeToElemData");
		return;
	}

	// count nr of points
	int NN = nodeMap->DataCount();

	// count nr of target points
	int NE = dom.Elements();
	int NP = 0;
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		NP += el.Nodes();
	}

	vector<double> srcData(NN);
	vector<double> trgData(NP);

	vector< vector<double> > mappedData(NP, vector<double>(9, 0.0));

	// loop over all the new nodes
	for (int l = 0; l < dataSize; ++l)
	{
		for (int i = 0; i < NN; ++i)
		{
			double vm = 0.0;
			switch (dataType)
			{
			case FEDataType::FE_DOUBLE: vm = nodeMap->value<double>(0, i); break;
			case FEDataType::FE_VEC3D:
			{
				vec3d v = nodeMap->value<vec3d>(0, i);
				if (l == 0) vm = v.x;
				if (l == 1) vm = v.y;
				if (l == 2) vm = v.z;
			}
			break;
			case FEDataType::FE_MAT3D:
			{
				int LUT[9][2] = { {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {1,2}, {2,0}, {2,1}, {2,2} };
				mat3d v = nodeMap->value<mat3d>(0, i);
				vm = v(LUT[l][0], LUT[l][1]);
			}
			break;
			case FEDataType::FE_MAT3DS:
			{
				int LUT[6][2] = { {0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2} };
				mat3ds v = nodeMap->value<mat3ds>(0, i);
				vm = v(LUT[l][0], LUT[l][1]);
			}
			break;
			default:
				assert(false);
			}
			srcData[i] = vm;
		}

		dataMapper->Map(trgData, [&srcData](int sourcePoint) {
			return srcData[sourcePoint];
		});

		for (int i = 0; i < NP; ++i)
		{
			mappedData[i][l] = trgData[i];
		}
	}

	// create new element set
	FEElementSet* eset = new FEElementSet(&fem);
	eset->Create(&dom);
	elemMap->Create(eset);

	// write mapped data to domain map
	int n = 0;
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);

		for (int k = 0; k < el.Nodes(); ++k)
		{
			vector<double>& vj = mappedData[n++];

			switch (dataType)
			{
			case FEDataType::FE_DOUBLE: elemMap->setValue(i, k, vj[0]); break;
			case FEDataType::FE_VEC3D:
			{
				vec3d v;
				v.x = vj[0];
				v.y = vj[1];
				v.z = vj[2];
				elemMap->setValue(i, k, v);
			}
			break;
			case FEDataType::FE_MAT3D:
			{
				mat3d v;
				v(0, 0) = vj[0]; v(0, 1) = vj[1]; v(0, 2) = vj[2];
				v(1, 0) = vj[3]; v(1, 1) = vj[4]; v(1, 2) = vj[5];
				v(2, 0) = vj[6]; v(2, 1) = vj[7]; v(2, 2) = vj[8];
				elemMap->setValue(i, k, v);
			}
			break;
			case FEDataType::FE_MAT3DS:
			{
				mat3ds v;
				v(0, 0) = vj[0];
				v(0, 1) = vj[1];
				v(0, 2) = vj[2];
				v(1, 1) = vj[3];
				v(1, 2) = vj[4];
				v(2, 2) = vj[5];
				elemMap->setValue(i, k, v);
			}
			break;
			default:
				assert(false);
			}
		}
	}
}

bool FEMMGRemesh::MMG::build_map_data(FEModel& fem)
{
	FEMesh& mesh = fem.GetMesh();

	m_nodeMapList.clear();
	m_nodeMapList.resize(mesh.Domains());

	// we'll need to store the original domain's nodal coordinates
	m_oldNodePos.resize(mesh.Domains());
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NN = dom.Nodes();

		vector<vec3d>& nodePos = m_oldNodePos[i];
		nodePos.resize(NN);
		for (int j = 0; j < NN; ++j) nodePos[j] = dom.Node(j).m_r0;
	}

	// loop over all domains
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);

		// create an element set (we need this for the domain map below)
		FEElementSet* elemSet = new FEElementSet(&fem);
		elemSet->Create(&dom);

		// write all material point data to a data stream
		DumpMemStream ar(fem);
		ar.Open(true, true);
		ar.WriteTypeInfo(true);

		// loop over all integration points
		int totalPoints = 0;
		for (int j = 0; j < dom.Elements(); ++j)
		{
			FEElement& el = dom.ElementRef(j);
			int nint = el.GaussPoints();
			for (int k = 0; k < nint; ++k)
			{
				FEMaterialPoint* mp = el.GetMaterialPoint(k);
				mp->Serialize(ar);
			}
			totalPoints += nint;
		}

		// figure out how much data was written for each material point
		size_t bytes = ar.bytesSerialized();

		size_t bytesPerPoint = bytes / totalPoints;
		assert((bytes%totalPoints) == 0);

		// re-open for reading
		ar.Open(false, true);

		// next, we need to figure out the datamaps for each data item
		vector<FEDomainMap*> mapList;
		DumpStream::DataBlock d;
		while (ar.bytesSerialized() < bytesPerPoint)
		{
			ar.readBlock(d);

			FEDomainMap* map = nullptr;
			switch (d.dataType())
			{
			case TypeID::TYPE_DOUBLE: map = new FEDomainMap(FEDataType::FE_DOUBLE, Storage_Fmt::FMT_MATPOINTS); break;
			case TypeID::TYPE_VEC3D: map = new FEDomainMap(FEDataType::FE_VEC3D, Storage_Fmt::FMT_MATPOINTS); break;
			case TypeID::TYPE_MAT3D: map = new FEDomainMap(FEDataType::FE_MAT3D, Storage_Fmt::FMT_MATPOINTS); break;
			case TypeID::TYPE_MAT3DS: map = new FEDomainMap(FEDataType::FE_MAT3DS, Storage_Fmt::FMT_MATPOINTS); break;
			default:
				assert(false);
				throw std::runtime_error("Error in mapping data.");
			}

			map->Create(elemSet, 0.0);

			mapList.push_back(map);
		}

		// rewind for processing
		ar.Open(false, true);

		for (int j = 0; j < dom.Elements(); ++j)
		{
			FEElement& el = dom.ElementRef(j);
			int nint = el.GaussPoints();
			for (int k = 0; k < nint; ++k)
			{
				int m = 0;
				size_t bytesRead = 0;
				while (bytesRead < bytesPerPoint)
				{
					size_t size0 = ar.bytesSerialized();
					bool b = ar.readBlock(d); assert(b);
					size_t size1 = ar.bytesSerialized();
					bytesRead += size1 - size0;

					FEDomainMap* map = mapList[m];

					switch (d.dataType())
					{
					case TypeID::TYPE_DOUBLE: { double v = d.value<double>(); map->setValue(j, k, v); } break;
					case TypeID::TYPE_MAT3D: { mat3d  v = d.value<mat3d >(); map->setValue(j, k, v); } break;
					case TypeID::TYPE_MAT3DS: { mat3ds v = d.value<mat3ds>(); map->setValue(j, k, v); } break;
					}

					m++;
					assert(m <= mapList.size());
				}
			}
		}

		// Now, we need to project all the data onto the nodes
		for (int j = 0; j < mapList.size(); ++j)
		{
			FEDomainMap* elemMap = mapList[j];
			FEDomainMap* nodeMap = createNodeDataMap(dom, elemMap);
			m_nodeMapList[i].push_back(nodeMap);
		}
	}

	// do the same thing for the mesh data
	m_meshDataList.clear();
	int dataMaps = mesh.DataMaps();
	for (int i = 0; i < dataMaps; ++i)
	{
		FEDomainMap* dmap = dynamic_cast<FEDomainMap*>(mesh.GetDataMap(i));
		if (dmap == nullptr) return false;
		if (dmap->DataType() != FEDataType::FE_DOUBLE) return false;

		const FEElementSet* elset = dmap->GetElementSet();
		const FEDomainList& domainList = elset->GetDomainList();
		if (domainList.Domains() != 1) return false;
		FEDomain& dom = const_cast<FEDomain&>(*domainList.GetDomain(0));

		// create a node data map of this domain map
		FEDomainMap* nodeMap = createNodeDataMap(dom, dmap);
		m_meshDataList.push_back(nodeMap);
	}

	return true;
}

void FEMMGRemesh::MMG::map_data(FEModel& fem)
{
	FEMesh& mesh = fem.GetMesh();

	// loop over all domains
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);

		// we need an element set for the domain maps below
		FEElementSet* elemSet = new FEElementSet(&fem);
		elemSet->Create(&dom);

		// build target node list
		vector<vec3d> trgPoints; trgPoints.reserve(dom.Elements());
		for (int i = 0; i < dom.Elements(); ++i)
		{
			FEElement& el = dom.ElementRef(i);
			int nint = el.GaussPoints();
			for (int j = 0; j < nint; ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				vec3d r = mp.m_r0;
				trgPoints.push_back(r);
			}
		}

		// set up mapper
		FEMeshDataInterpolator* mapper = nullptr;
		switch (m_mmgRemesh->m_transferMethod)
		{
		case TRANSFER_SHAPE:
		{
			FEDomain* oldDomain = &m_meshCopy->Domain(i);
			FEDomainShapeInterpolator* dsm = new FEDomainShapeInterpolator(oldDomain);
			dsm->SetTargetPoints(trgPoints);
			mapper = dsm;
		}
		break;
		case TRANSFER_MLQ:
		{
			FELeastSquaresInterpolator* MLQ = new FELeastSquaresInterpolator;
			MLQ->SetNearestNeighborCount(m_mmgRemesh->m_nnc);
			MLQ->SetSourcePoints(m_oldNodePos[i]);
			MLQ->SetTargetPoints(trgPoints);
			mapper = MLQ;
		}
		break;
		default:
			assert(false);
			return;
		}
		if (mapper->Init() == false)
		{
			assert(false);
			throw std::runtime_error("Failed to initialize LLQ");
		}

		// loop over all the domain maps
		std::vector<FEDomainMap*>& nodeMap_i = m_nodeMapList[i];
		int mapCount = nodeMap_i.size();
		vector<FEDomainMap*> elemMapList(mapCount);
		for (int j = 0; j < mapCount; ++j)
		{
			FEDomainMap* nodeMap = nodeMap_i[j];

			// map node data to integration points
			FEDomainMap* elemMap = createElemDataMap(fem, dom, m_oldNodePos[i], nodeMap, mapper);

			elemMapList[j] = elemMap;
		}

		// now we need to reconstruct the data stream
		DumpMemStream ar(fem);
		ar.Open(true, true);

		for (int j = 0; j < dom.Elements(); ++j)
		{
			FEElement& el = dom.ElementRef(j);
			int nint = el.GaussPoints();

			for (int k = 0; k < nint; ++k)
			{
				for (int l = 0; l < mapCount; ++l)
				{
					FEDomainMap* map = elemMapList[l];
					switch (map->DataType())
					{
					case FEDataType::FE_DOUBLE: { double v = map->value<double>(j, k); ar << v; } break;
					case FEDataType::FE_VEC3D : { vec3d  v = map->value<vec3d >(j, k); ar << v; } break;
					case FEDataType::FE_MAT3D : { mat3d  v = map->value<mat3d >(j, k); ar << v; } break;
					case FEDataType::FE_MAT3DS: { mat3ds v = map->value<mat3ds>(j, k); ar << v; } break;
					default:
						assert(false);
					}
				}
			}
		}

		// time to serialize everything back to the new integration points
		ar.Open(false, true);
		for (int j = 0; j < dom.Elements(); ++j)
		{
			FEElement& el = dom.ElementRef(j);
			int nint = el.GaussPoints();

			for (int k = 0; k < nint; ++k)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(k);
				mp.Serialize(ar);
			}
		}
	}

	// map mesh data
	for (int i = 0; i < m_meshDataList.size(); ++i)
	{
		FEDomainMap* elemMap = dynamic_cast<FEDomainMap*>(mesh.GetDataMap(i));
		FEDomainMap* nodeMap = m_meshDataList[i];

		// build target points.
		const FEDomainList& domainList = elemMap->GetElementSet()->GetDomainList();
		FEDomain& dom = const_cast<FEDomain&>(*domainList.GetDomain(0));
		int NE = dom.Elements();
		vector<vec3d> trgPoints; trgPoints.reserve(NE);
		for (int n = 0; n < NE; ++n)
		{
			FEElement& el = dom.ElementRef(n);
			int ne = el.Nodes();
			for (int l = 0; l < ne; ++l)
			{
				vec3d r = mesh.Node(el.m_node[l]).m_r0;
				trgPoints.push_back(r);
			}
		}

		// set up mapper
		FEMeshDataInterpolator* mapper = nullptr;
		switch (m_mmgRemesh->m_transferMethod)
		{
		case TRANSFER_SHAPE:
		{
			FEDomain* oldDomain = m_meshCopy->FindDomain(dom.GetName());
			FEDomainShapeInterpolator* dsm = new FEDomainShapeInterpolator(oldDomain);
			dsm->SetTargetPoints(trgPoints);
			mapper = dsm;
		}
		break;
		case TRANSFER_MLQ:
		{
			FELeastSquaresInterpolator* MLQ = new FELeastSquaresInterpolator;
			MLQ->SetNearestNeighborCount(m_mmgRemesh->m_nnc);
			MLQ->SetSourcePoints(m_oldNodePos[i]);
			MLQ->SetTargetPoints(trgPoints);
			mapper = MLQ;
		}
		break;
		default:
			assert(false);
			return;
		}
		if (mapper->Init() == false)
		{
			assert(false);
			throw std::runtime_error("Failed to initialize LLQ");
		}

		// map node data to integration points
		NodeToElemData(fem, dom, nodeMap, elemMap, mapper);
	}
}

#endif
