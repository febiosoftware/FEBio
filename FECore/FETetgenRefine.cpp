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
#include "FETetgenRefine.h"
#include "FEModel.h"
#include "FEMeshTopo.h"
#include "FESolidDomain.h"
#include "FESurface.h"
#include "FESurfaceLoad.h"
#include "FEBoundaryCondition.h"
#include "FESurfacePairConstraint.h"
#include "log.h"

#ifdef TETLIBRARY
#include <tetgen.h>
#endif

struct TETGENOPTIONS
{
	double	h;
	double	q;
};

bool build_tetgen_remesh(FEMeshTopo& topo, tetgenio& in, vector<int>& elemSelection, TETGENOPTIONS& ops);

BEGIN_FECORE_CLASS(FETetgenRefine, FERefineMesh)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_q, "q");
	ADD_PARAMETER(m_tol, "tol");
	ADD_PARAMETER(m_splitFaces, "split_faces");
	ADD_PARAMETER(m_maxiter, "max_iters");
	ADD_PARAMETER(m_maxelem, "max_elems");
	ADD_PARAMETER(m_resetMesh, "reset_mesh");
	ADD_PROPERTY(m_criterion, "criterion", 0);
END_FECORE_CLASS();

FETetgenRefine::FETetgenRefine(FEModel* fem) : FERefineMesh(fem)
{
	m_scale = 0.5;
	m_q = 0.0;
	m_tol = 0.0;
	m_splitFaces = true;
	m_maxiter = 1;
	m_maxelem = 0;
	m_criterion = nullptr;
	m_resetMesh = false;
}

bool FETetgenRefine::Apply(int iteration)
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
	if (DoTetRefinement(fem) == false)
	{
		feLogError("Nothing to do.");
		return true;
	}

	// reactivate the model
	// reactivate BCs
	for (int i = 0; i < fem.BoundaryConditions(); ++i)
	{
		FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
		if (bc.IsActive()) bc.Activate();
	}

	// update surface loads 
	for (int i = 0; i < fem.SurfaceLoads(); ++i)
	{
		FESurfaceLoad& sl = *fem.SurfaceLoad(i);
		sl.SetSurface(&sl.GetSurface());
		if (sl.IsActive()) sl.Activate();
	}
	// update surface interactions
	for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
	{
		FESurfacePairConstraint& ci = *fem.SurfacePairConstraint(i);
		if (ci.IsActive()) ci.Activate();
	}

	// all done
	return false;
}

FESolidElement* FindElement(FEMesh& mesh, const vec3d& y, double r[3])
{
	FESolidElement* el = nullptr;
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FESolidDomain* dom = dynamic_cast<FESolidDomain*>(&mesh.Domain(i));
		el = dom->FindReferenceElement(y, r);
		if (el) return el;
	}
	return nullptr;
}

bool FETetgenRefine::DoTetRefinement(FEModel& fem)
{
#ifdef TETLIBRARY

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// make sure it is a tet4 mesh
	if (mesh.IsType(ET_TET4)==false)
	{
		feLogError("This is not a TET4 mesh");
		return 0;
	}

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

	// allocate tetgen structures
	tetgenio in, out;
	in.initialize();
	out.initialize();

	FEMeshTopo& topo = *m_topo;

	// build the tetgen structure
	TETGENOPTIONS tgops;
	tgops.h = m_scale;
	if (build_tetgen_remesh(*m_topo, in, m_elemList, tgops) == false) return 0;

	// set the parameters
	double q = m_q;
	bool bsplit = m_splitFaces;

	char sz[64] = { 0 };
	char* ch = sz;
	sprintf(ch, "r");  ++ch;
	if (q   > 0) sprintf(ch, "q%lg", q); ch += strlen(ch);
	sprintf(ch, "a"); ch += strlen(ch);
	if (m_tol > 0) sprintf(ch, "T%lg", m_tol); ch += strlen(ch);
	if (bsplit == false) sprintf(ch, "Y"); ++ch;

	// create a tet mesh
	tetrahedralize(sz, &in, &out);
	assert(out.numberofcorners == 4);

	int nodes = out.numberofpoints;
	int elems = out.numberoftetrahedra;
	int faces = out.numberoftrifaces;

	// reallocate nodes
	int N0 = mesh.Nodes();
	int N1 = nodes;
	mesh.AddNodes(N1 - N0);

	// copy nodes
	for (int i = N0; i<nodes; ++i)
	{
		vec3d r;
		r.x = out.pointlist[3*i    ];
		r.y = out.pointlist[3*i + 1];
		r.z = out.pointlist[3*i + 2];

		FENode& node = mesh.Node(i);
		node.m_r0 = r;
	}

	// assign dofs to new nodes
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	for (int i = N0; i < nodes; ++i)
	{
		FENode& node = mesh.Node(i);
		node.SetDOFS(MAX_DOFS);
		node.m_rt = node.m_r0;
		for (int j = 0; j < node.m_ID.size(); ++j) {
			node.set(j, 0.0);
		}
	}

	if (m_resetMesh)
	{
		for (int i = 0; i < N0; ++i)
		{
			FENode& node = mesh.Node(i);
			node.m_rt = node.m_r0;
			for (int j = 0; j < node.m_ID.size(); ++j) {
				node.set(j, 0.0);
			}
		}
	}
	else
	{
		// update solution
		for (int i = N0; i < nodes; ++i)
		{
			FENode& nodei = mesh.Node(i);

			// find the element in which this node is
			double r[3] = { 0 };
			FESolidElement* el = FindElement(mesh, nodei.m_r0, r);
			if (el == nullptr)
			{
				assert(false);
				return false;
			}

			// get the nodal coordinates
			vec3d rt[FEElement::MAX_NODES];
			for (int j = 0; j < el->Nodes(); ++j) rt[j] = mesh.Node(el->m_node[j]).m_rt;

			nodei.m_rt = el->evaluate(rt, r[0], r[1], r[2]);

			// update values
			for (int l = 0; l < MAX_DOFS; ++l)
			{
				double v[FEElement::MAX_NODES] = { 0 };
				for (int j = 0; j < el->Nodes(); ++j) v[j] = mesh.Node(el->m_node[j]).get(l);
				double vl = el->evaluate(v, r[0], r[1], r[2]);
				nodei.set(l, vl);
			}
			nodei.UpdateValues();
		}
	}

	// count the number of edges with a marker > 1
	int edges = 0;
	for (int i=0; i<out.numberofedges; ++i)
	{
		if (out.edgemarkerlist[i] > 1) edges++;
	}

	// recreate domains
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(i));

		int nelems = 0;
		for (int j = 0; j < elems; ++j)
		{
			if (out.tetrahedronattributelist[j] == i) nelems++;
		}

		dom.Create(nelems, FE_TET4G1);
		int n = 0;
		for (int j = 0; j < elems; ++j)
		{
			if (out.tetrahedronattributelist[j] == i)
			{
				FESolidElement& el = dom.Element(n++);
				el.m_node[0] = out.tetrahedronlist[4 * j];
				el.m_node[1] = out.tetrahedronlist[4 * j + 1];
				el.m_node[2] = out.tetrahedronlist[4 * j + 2];
				el.m_node[3] = out.tetrahedronlist[4 * j + 3];
			}
		}
	}

/*	// copy faces
	for (int i=0; i<faces; ++i)
	{
		FEFace& f = pmesh->Face(i);
		f.SetType(FE_FACE_TRI3);
		f.n[0] = out.trifacelist[3 * i + 2];
		f.n[1] = out.trifacelist[3*i+1];
		f.n[2] = out.trifacelist[3*i  ];
		f.n[3] = f.n[2];
		f.m_gid  = gid[out.trifacemarkerlist[i]];
		f.m_sid  = sid[out.trifacemarkerlist[i]];
	}
*/
	// re-init domains
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.CreateMaterialPointData();
		dom.Init();
		dom.Activate();
	}

	vector<int> fm(out.numberoftrifaces);
	for (int i = 0; i < out.numberoftrifaces; ++i) fm[i] = out.trifacemarkerlist[i];

	// recreate surfaces
	int faceMark = 1;
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		FESurface& surf = mesh.Surface(i);

		// count faces
		int nfaces = 0;
		for (int j = 0; j < out.numberoftrifaces; ++j)
		{
			if (out.trifacemarkerlist[j] == faceMark) nfaces++;
		}
		assert(nfaces > 0);
		surf.Create(nfaces);
		int n = 0;
		for (int j = 0; j < out.numberoftrifaces; ++j)
		{
			if (out.trifacemarkerlist[j] == faceMark)
			{
				FESurfaceElement& face = surf.Element(n++);
				face.SetType(FE_TRI3G3);
				face.m_node[0] = out.trifacelist[3*j+2];
				face.m_node[1] = out.trifacelist[3*j+1];
				face.m_node[2] = out.trifacelist[3*j  ];
			}
		}
		surf.Init();

		faceMark++;
	}

	// update nodesets
	for (int i = 0; i < mesh.NodeSets(); ++i)
	{
		FENodeSet& nset = *mesh.NodeSet(i);
		vector<int> nodeTags(mesh.Nodes(), 0);
		for (int j = 0; j < out.numberoftrifaces; ++j)
		{
			if (out.trifacemarkerlist[j] == faceMark)
			{
				int* fn = out.trifacelist + 3 * j;
				nodeTags[fn[0]] = 1;
				nodeTags[fn[1]] = 1;
				nodeTags[fn[2]] = 1;
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

	return true;
#else 
	return false;
#endif // TETLIBRARY
}

//-----------------------------------------------------------------------------
bool build_tetgen_remesh(FEMeshTopo& topo, tetgenio& in, vector<int>& elemSelection, TETGENOPTIONS& ops )
{
#ifdef TETLIBRARY
	FEMesh& mesh = *topo.GetMesh();

	// make sure this is a tetmesh
	assert(mesh.IsType(ET_TET4));

	// all indices start from 0
	in.firstnumber = 0;

	// allocate nodes
	int nodes = mesh.Nodes();
	in.numberofpoints = nodes;
	in.pointlist = new REAL[3 * nodes];
	for (int i = 0; i<nodes; ++i)
	{
		vec3d& r = mesh.Node(i).m_r0;
		in.pointlist[3*i    ] = r.x;
		in.pointlist[3*i + 1] = r.y;
		in.pointlist[3*i + 2] = r.z;
	}

	// build the element list
	int elems = topo.Elements();
	in.numberoftetrahedra = elems;
	in.tetrahedronlist = new int[elems * 4];
	int n = 0;
	for (int i = 0; i<mesh.Domains(); ++i)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(i));
		for (int j = 0; j < dom.Elements(); ++j, ++n)
		{
			FESolidElement& el = dom.Element(j);
			in.tetrahedronlist[4 * n    ] = el.m_node[0];
			in.tetrahedronlist[4 * n + 1] = el.m_node[1];
			in.tetrahedronlist[4 * n + 2] = el.m_node[2];
			in.tetrahedronlist[4 * n + 3] = el.m_node[3];
		}
	}

	// set the element attributes
	in.numberoftetrahedronattributes = 1;
	in.tetrahedronattributelist = new double[elems];
	n = 0;
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(i));
		for (int j = 0; j < dom.Elements(); ++j, ++n)
		{
			in.tetrahedronattributelist[n] = i;
		}
	}

	// build the facet list
	int faces = topo.SurfaceFaces();
	in.numberoftrifaces = faces;
	in.trifacelist = new int[faces * 3];
	for (int i = 0, n = 0; i<faces; ++i)
	{
		const FEFaceList::FACE& f = topo.SurfaceFace(i);
		assert(f.ntype == 3);
		in.trifacelist[3*i    ] = f.node[0];
		in.trifacelist[3*i + 1] = f.node[1];
		in.trifacelist[3*i + 2] = f.node[2];
	}

	// set the facet markers
	in.trifacemarkerlist = new int[faces];
	for (int i = 0; i<faces; ++i) in.trifacemarkerlist[i] = 0;

	int faceMark = 1;
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		FESurface& surf = mesh.Surface(i);
		vector<int> faceIndexList = topo.SurfaceFaceIndexList(surf);
		for (int j = 0; j < faceIndexList.size(); ++j)
		{
			in.trifacemarkerlist[faceIndexList[j]] = faceMark;
		}
		faceMark++;
	}

	// for node sets we are going to create artificial surfaces
	// (skip nodeset 0 since that is all the nodes)
	vector<int> nodeTags;
	for (int i = 0; i < mesh.NodeSets(); ++i)
	{
		FENodeSet& nset = *mesh.NodeSet(i);
		nodeTags.assign(mesh.Nodes(), 0);
		for (int j = 0; j < nset.Size(); ++j) nodeTags[nset[j]] = 2;

		// see if this is indeed a surface node set
		for (int j = 0; j < in.numberoftrifaces; ++j)
		{
			int* fn = in.trifacelist + 3 * j;
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
			for (int j = 0; j < in.numberoftrifaces; ++j)
			{
				int* fn = in.trifacelist + 3 * j;
				if ((nodeTags[fn[0]] == 1) && (nodeTags[fn[1]] == 1) && (nodeTags[fn[2]] == 1))
				{
					assert(in.trifacemarkerlist[j] == 0);
					in.trifacemarkerlist[j] = faceMark;
				}
			}
		}
		faceMark++;
	}

	double h = ops.h;

	// build the facet constraint list
/*	in.numberoffacetconstraints = 0;
	in.facetconstraintlist = new double[2 * 1];
	double a = h*h;
	for (int i = 0; i<1; ++i)
	{
		const FEFaceList::FACE& f = topo.Face(i);
		in.facetconstraintlist[2 * i] = 0;
		in.facetconstraintlist[2 * i + 1] = a;
	}
*/
/*
	int nfeather = GetIntValue(3);
	if (nfeather > 0)
	{
		vector<double> area; area.assign(faces, 0.0);
		for (int i = 0; i<faces; ++i)
		{
			FEFace& f = pm->Face(i);
			area[i] = FEMeshMetrics::SurfaceArea(*pm, f);
		}

		for (int i = 0; i<faces; ++i)
		{
			FEFace& f = pm->Face(i);
			FEElement& el = pm->Element(f.m_elem[0]);
			if (el.IsSelected() || f.IsSelected()) f.m_ntag = 1; else f.m_ntag = 0;
		}

		for (int n = 0; n<nfeather; ++n)
		{
			double w = (n + 1.0) / (nfeather + 1.0);
			w *= w;
			for (int i = 0; i<faces; ++i)
			{
				FEFace& f = pm->Face(i);
				if (f.m_ntag == 1)
				{
					int nf = 3; assert(f.Nodes() == 3);
					for (int j = 0; j<nf; ++j)
					{
						FEFace* f2 = pm->FacePtr(f.m_nbr[j]);
						if (f2 && (f2->m_ntag == 0))
						{
							in.facetconstraintlist[2 * f.m_nbr[j] + 1] = a*(1.0 - w) + w*area[f.m_nbr[j]];
							f2->m_ntag = 2;
						}
					}
				}
			}
			for (int i = 0; i<faces; ++i)
			{
				FEFace& f = pm->Face(i);
				if (f.m_ntag == 2) f.m_ntag = 1;
			}
		}
	}
*/

	// build the element volume list
	if (h > 0)
	{
		in.tetrahedronvolumelist = new REAL[elems];
		for (int i = 0; i<elems; ++i)
		{
			if (elemSelection[i] == 1)
			{
				// calculate volume of tet
				double V = mesh.ElementVolume(*topo.Element(i));
				in.tetrahedronvolumelist[i] = V*h;
			}
			else
				in.tetrahedronvolumelist[i] = 0;
		}

/*		if (nfeather > 0)
		{
			vector<double> evol; evol.assign(elems, 0.0);
			for (int i = 0; i<elems; ++i)
			{
				FEElement& el = pm->Element(i);
				evol[i] = FEMeshMetrics::ElementVolume(*pm, el);
			}

			for (int i = 0; i<elems; ++i)
			{
				FEElement& el = pm->Element(i);
				if (el.IsSelected()) el.m_ntag = 1; else el.m_ntag = 0;
			}

			for (int n = 0; n<nfeather; ++n)
			{
				double w = (n + 1.0) / (nfeather + 1.0);
				w *= w;
				for (int i = 0; i<elems; ++i)
				{
					FEElement& el = pm->Element(i);
					if (el.m_ntag == 1)
					{
						int nf = el.Faces();
						for (int j = 0; j<nf; ++j)
						{
							FEElement* pe2 = pm->ElementPtr(el.m_nbr[j]);
							if (pe2 && (pe2->m_ntag == 0))
							{
								in.tetrahedronvolumelist[el.m_nbr[j]] = a*(1.0 - w) + w*evol[el.m_nbr[j]];
								pe2->m_ntag = 2;
							}
						}
					}
				}
				for (int i = 0; i<elems; ++i)
				{
					FEElement& el = pm->Element(i);
					if (el.m_ntag == 2) el.m_ntag = 1;
				}
			}
		}
		*/
	}

	// define the edges
/*	in.numberofedges = pm->Edges();
	in.edgelist = new int[2 * in.numberofedges];
	in.edgemarkerlist = new int[in.numberofedges];
	for (int i = 0; i<in.numberofedges; ++i)
	{
		FEEdge& e = pm->Edge(i);
		in.edgelist[2 * i] = e.n[0];
		in.edgelist[2 * i + 1] = e.n[1];
		in.edgemarkerlist[i] = e.m_gid + 2;
	}
*/
	return true;

#else
	return false;
#endif
}
