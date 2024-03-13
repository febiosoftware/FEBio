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
#include "FEParabolicMap.h"
#include "FESurfaceMap.h"
#include "FESurface.h"
#include "FEElemElemList.h"
#include "log.h"
#include "SparseMatrix.h"
#include "LinearSolver.h"
#include "FEGlobalMatrix.h"
#include "FEModel.h"

BEGIN_FECORE_CLASS(FEParabolicMap, FEFaceDataGenerator)
	ADD_PARAMETER(m_scale, "value");
END_FECORE_CLASS();

FEParabolicMap::FEParabolicMap(FEModel* fem) : FEFaceDataGenerator(fem), m_dofs(fem)
{
	m_scale = 1.0;
}

FEParabolicMap::~FEParabolicMap()
{

}

void FEParabolicMap::SetDOFConstraint(const FEDofList& dofs)
{
	m_dofs = dofs;
}

FEDataMap* FEParabolicMap::Generate()
{
	const FEFacetSet& facetSet = *GetFacetSet();
	FESurfaceMap* map = new FESurfaceMap(FEDataType::FE_DOUBLE);
	map->Create(&facetSet, 0.0, FMT_NODE);

	// create temporary surface
	FESurface surf(GetFEModel());
	surf.Create(facetSet);
	surf.InitSurface();

	// find surface boundary nodes
	FEElemElemList EEL;
	EEL.Create(&surf);

	vector<bool> boundary(surf.Nodes(), false);
	for (int i = 0; i<surf.Elements(); ++i) {
		FESurfaceElement& el = surf.Element(i);
		for (int j = 0; j<el.facet_edges(); ++j) {
			FEElement* nel = EEL.Neighbor(i, j);
			if (nel == nullptr) {
				int en[3] = { -1,-1,-1 };
				el.facet_edge(j, en);
				boundary[en[0]] = true;
				boundary[en[1]] = true;
				if (en[2] > -1) boundary[en[2]] = true;
			}
		}
	}

	// Apply dof constraints
	if (m_dofs.IsEmpty() == false)
	{
		// only consider nodes with fixed dofs as boundary nodes
		for (int i = 0; i < surf.Nodes(); ++i)
			if (boundary[i]) {
				FENode& node = surf.Node(i);

				bool b = false;
				for (int j = 0; j < m_dofs.Size(); ++j)
				{
					if (node.get_bc(m_dofs[j]) != DOF_FIXED) b = true;
				}

				if (b) boundary[i] = false;
			}
	}

	// count number of non-boundary nodes
	int neq = 0;
	vector<int> glm(surf.Nodes(), -1);
	for (int i = 0; i< surf.Nodes(); ++i)
		if (!boundary[i]) glm[i] = neq++;
	if (neq == 0)
	{
		feLogError("Unable to set parabolic map\n");
		return nullptr;
	}

	// create a linear solver
	LinearSolver*		plinsolve;	//!< the linear solver
	FEGlobalMatrix*		pK;			//!< stiffness matrix
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	plinsolve = fecore_new<LinearSolver>("skyline", nullptr);
	if (plinsolve == 0)
	{
		feLogError("Unknown solver type selected\n");
		return nullptr;
	}

	SparseMatrix* pS = plinsolve->CreateSparseMatrix(REAL_SYMMETRIC);
	pK = new FEGlobalMatrix(pS);
	if (pK == 0)
	{
		feLogError("Failed allocating stiffness matrix\n\n");
		return nullptr;
	}
	// build matrix profile for normal velocity at non-boundary nodes
	pK->build_begin(neq);
	for (int i = 0; i< surf.Elements(); ++i) {
		FESurfaceElement& el = surf.Element(i);
		vector<int> elm(el.Nodes(), -1);
		for (int j = 0; j<el.Nodes(); ++j)
			elm[j] = glm[el.m_lnode[j]];
		pK->build_add(elm);
	}
	pK->build_end();
	pS->Zero();

	// create global vector
	vector<double> v;           //!< solution
	vector<double> rhs;         //!< right-hand-side
	vector<double> Fr;          //!< reaction forces
	v.assign(neq, 0);
	rhs.assign(neq, 0);
	Fr.assign(neq, 0);
	FEModel pfem;
	FEGlobalVector pR(pfem, rhs, Fr);

	// calculate the global matrix and vector
	FEElementMatrix ke;
	vector<double> fe;
	vector<int> lm;

	for (int m = 0; m< surf.Elements(); ++m)
	{
		// get the surface element
		FESurfaceElement& el = surf.Element(m);

		int neln = el.Nodes();

		// get the element stiffness matrix
		ke.resize(neln, neln);
		lm.resize(neln);
		fe.resize(neln);
		vector<vec3d> gradN(neln);

		// calculate stiffness
		int nint = el.GaussPoints();

		// gauss weights
		double* w = el.GaussWeights();

		// nodal coordinates
		FEMesh& mesh = *surf.GetMesh();
		vec3d rt[FEElement::MAX_NODES];
		for (int j = 0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;

		// repeat over integration points
		ke.zero();
		zero(fe);
		vec3d gcnt[2];
		for (int n = 0; n<nint; ++n)
		{
			double* N = el.H(n);
			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);
			surf.ContraBaseVectors(el, n, gcnt);

			vec3d dxr(0, 0, 0), dxs(0, 0, 0);
			for (int i = 0; i<neln; ++i)
			{
				dxr += rt[i] * Gr[i];
				dxs += rt[i] * Gs[i];
				gradN[i] = gcnt[0] * Gr[i] + gcnt[1] * Gs[i];
			}

			double da = (dxr ^ dxs).norm();

			// calculate stiffness component
			for (int i = 0; i<neln; ++i) {
				fe[i] += N[i] * w[n] * da;
				for (int j = 0; j<neln; ++j)
					ke[i][j] += (gradN[i] * gradN[j])*w[n] * da;
			}
		}

		// get the element's LM vector
		for (int j = 0; j<el.Nodes(); ++j)
			lm[j] = glm[el.m_lnode[j]];

		// assemble element matrix in global stiffness matrix
		ke.SetIndices(lm);
		pK->Assemble(ke);
		pR.Assemble(lm, fe);
	}

	// solve linear system
	plinsolve->PreProcess();
	plinsolve->Factor();
	if (plinsolve->BackSolve(v, rhs) == false)
	{
		feLogError("Unable to solve for parabolic field\n");
		return nullptr;
	}
	plinsolve->Destroy();

	// set the nodal normal velocity scale factors
	vector<double> VN(surf.Nodes(), 0.0);
	for (int i = 0; i< surf.Nodes(); ++i) {
		if (glm[i] == -1) VN[i] = 0;
		else VN[i] = v[glm[i]];
	}

	// evaluate net area and volumetric flow rate
	double A = 0, Q = 0;
	for (int m = 0; m< surf.Elements(); ++m)
	{
		// get the surface element
		FESurfaceElement& el = surf.Element(m);

		int neln = el.Nodes();
		int nint = el.GaussPoints();
		double* w = el.GaussWeights();

		// nodal coordinates
		FEMesh& mesh = *surf.GetMesh();
		vec3d rt[FEElement::MAX_NODES];
		for (int j = 0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;

		// repeat over integration points
		for (int n = 0; n<nint; ++n)
		{
			double* N = el.H(n);
			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);

			double vn = 0;
			vec3d dxr(0, 0, 0), dxs(0, 0, 0);
			for (int i = 0; i<neln; ++i)
			{
				vn += N[i] * VN[el.m_lnode[i]];
				dxr += rt[i] * Gr[i];
				dxs += rt[i] * Gs[i];
			}

			double da = (dxr ^ dxs).norm();

			for (int i = 0; i<neln; ++i) {
				A += N[i] * w[n] * da;
				Q += N[i] * vn*w[n] * da;
			}
		}
	}

	// normalize nodal velocity cards
	double vbar = Q / A;
	for (int i = 0; i< surf.Nodes(); ++i) VN[i] /= vbar;

	// assign nodal values to surface map
	map->set<double>(0.0);
	for (int i = 0; i < surf.Nodes(); ++i)
	{
		map->set<double>(i, VN[i] * m_scale);
	}

	return map;
}
