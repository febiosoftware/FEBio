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
#include "log.h"

#ifdef TETLIBRARY
#include <tetgen.h>
#endif

struct TETGENOPTIONS
{
	double	h;
	double	q;
};

bool build_tetgen_remesh(FEMeshTopo& topo, tetgenio& in, TETGENOPTIONS& ops);

BEGIN_FECORE_CLASS(FETetgenRefine, FERefineMesh)
	ADD_PARAMETER(m_h, "h");
	ADD_PARAMETER(m_q, "q");
	ADD_PARAMETER(m_tol, "tol");
	ADD_PARAMETER(m_splitFaces, "split_faces");
	ADD_PARAMETER(m_maxiter, "max_iters");
	ADD_PARAMETER(m_maxelem, "max_elems");
	ADD_PROPERTY(m_criterion, "criterion", 0);
END_FECORE_CLASS();

FETetgenRefine::FETetgenRefine(FEModel* fem) : FERefineMesh(fem)
{
	m_h = 0.0;
	m_q = 0.0;
	m_tol = 0.0;
	m_splitFaces = true;
	m_maxiter = 1;
	m_maxelem = 0;
	m_criterion = nullptr;
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
		feLogError("Something went wrong applying tetgen adaptor.");
		return true;
	}

	// all done
	return (iteration > 0);
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

	// allocate tetgen structures
	tetgenio in, out;
	in.initialize();
	out.initialize();

	FEMeshTopo& topo = *m_topo;
	int NF = topo.Faces();
	vector<int> gid(NF), sid(NF);
	for (int i=0; i<NF; ++i) 
	{
		const FEFaceList::FACE& face = topo.Face(i);
//		gid[i] = face.m_gid;
	}

	// build the tetgen structure
	TETGENOPTIONS tgops;
	tgops.h = m_h;
	if (build_tetgen_remesh(*m_topo, in, tgops) == false) return 0;

	// convert element size to volume
	double h = m_h;
	double a = h*h*h;

	// set the parameters
	double q = m_q;
	bool bsplit = m_splitFaces;

	char sz[64] = { 0 };
	char* ch = sz;
	sprintf(ch, "r");  ++ch;
	if (q   > 0) sprintf(ch, "q%lg", q); ch += strlen(ch);
	if (h   > 0) sprintf(ch, "a"); ch += strlen(ch);
	if (m_tol > 0) sprintf(ch, "T%lg", m_tol); ch += strlen(ch);
	if (bsplit == false) sprintf(ch, "Y"); ++ch;

	// create a tet mesh
	tetrahedralize(sz, &in, &out);
	assert(out.numberofcorners == 4);

	int nodes = out.numberofpoints;
	int elems = out.numberoftetrahedra;
	int faces = out.numberoftrifaces;

	// reallocate nodes
	mesh.CreateNodes(nodes);

	// copy nodes
	for (int i = 0; i<nodes; ++i)
	{
		vec3d r;
		r.x = out.pointlist[3*i    ];
		r.y = out.pointlist[3*i + 1];
		r.z = out.pointlist[3*i + 2];

		FENode& node = mesh.Node(i);
		node.m_r0 = node.m_rt = node.m_rp = r;
	}

	// assign dofs to new nodes
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	mesh.SetDOFS(MAX_DOFS);
	for (int i = 0; i < nodes; ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j = 0; j < node.m_ID.size(); ++j) {
			node.set(j, 0.0);
		}
		node.UpdateValues();
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
		dom.Create(elems, FE_TET4G1);
		for (int j = 0; j < elems; ++j)
		{
			FESolidElement& el = dom.Element(j);
			el.m_node[0] = out.tetrahedronlist[4*j    ];
			el.m_node[1] = out.tetrahedronlist[4*j + 1];
			el.m_node[2] = out.tetrahedronlist[4*j + 2];
			el.m_node[3] = out.tetrahedronlist[4*j + 3];
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

	return true;
#else 
	return false;
#endif // TETLIBRARY
}

//-----------------------------------------------------------------------------
bool build_tetgen_remesh(FEMeshTopo& topo, tetgenio& in, TETGENOPTIONS& ops)
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
	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(0));
	int elems = dom.Elements();
	in.numberoftetrahedra = elems;
	in.tetrahedronlist = new int[elems * 4];
	for (int i = 0; i<elems; ++i)
	{
		FESolidElement& el = dom.Element(i);
		in.tetrahedronlist[4*i    ] = el.m_node[0];
		in.tetrahedronlist[4*i + 1] = el.m_node[1];
		in.tetrahedronlist[4*i + 2] = el.m_node[2];
		in.tetrahedronlist[4*i + 3] = el.m_node[3];
	}

	// build the facet list
	int faces = topo.Faces();
	in.numberoftrifaces = faces;
	in.trifacelist = new int[faces * 3];
	for (int i = 0, n = 0; i<faces; ++i)
	{
		const FEFaceList::FACE& f = topo.Face(i);
		assert(f.ntype == 3);
		in.trifacelist[3*i    ] = f.node[0];
		in.trifacelist[3*i + 1] = f.node[1];
		in.trifacelist[3*i + 2] = f.node[2];
	}

	// set the facet markers
	in.trifacemarkerlist = new int[faces];
	for (int i = 0; i<faces; ++i) in.trifacemarkerlist[i] = 0;

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
		double a = h*h*h;
		in.tetrahedronvolumelist = new REAL[elems];
		for (int i = 0; i<elems; ++i)
		{
			in.tetrahedronvolumelist[i] = a;
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
