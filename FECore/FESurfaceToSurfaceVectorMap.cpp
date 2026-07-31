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
#include "FESurfaceToSurfaceVectorMap.h"
#include "FEMesh.h"
#include "FESurface.h"
#include "FEDomainMap.h"
#include "FECoreKernel.h"
#include "LinearSolver.h"
#include "FEGlobalMatrix.h"
#include "FELinearSystem.h"
#include "FESolidDomain.h"
#include "FEModel.h"
#include "log.h"

BEGIN_FECORE_CLASS(FESurfaceToSurfaceVectorMap, FEElemDataGenerator)
	ADD_PROPERTY(m_surf[0], "surface1", FEProperty::Reference);
	ADD_PROPERTY(m_surf[1], "surface2", FEProperty::Reference);
	ADD_PARAMETER(m_normal, "normal");
	ADD_PARAMETER(m_cross, "cross");
	ADD_PARAMETER(m_inAngle, "in_angle");
	ADD_PARAMETER(m_outAngle, "out_angle");
	ADD_PARAMETER(m_smoothIters, "smooth_iters");
END_FECORE_CLASS();

FESurfaceToSurfaceVectorMap::FESurfaceToSurfaceVectorMap(FEModel* fem) : FEElemDataGenerator(fem)
{
	m_surf[0] = nullptr;
	m_surf[1] = nullptr;
	m_normal = vec3d(0, 0, 1);

	m_cross = false;
	m_inAngle = 0;
	m_outAngle = 0;
}

FESurfaceToSurfaceVectorMap::~FESurfaceToSurfaceVectorMap()
{
}

bool FESurfaceToSurfaceVectorMap::Init()
{
	if ((m_surf[0] == nullptr) || (m_surf[1] == nullptr)) return false;
	return FEMeshDataGenerator::Init();
}

FEDataMap* FESurfaceToSurfaceVectorMap::Generate()
{
	FEElementSet* elset = GetElementSet();
	if (elset == nullptr) return nullptr;

	FEMesh& mesh = GetMesh();

	int NN = mesh.Nodes();
	vector<int> bn_mesh(NN, 0);
	vector<double> val_mesh(NN, 0.0);

	for (int i = 0; i < 2; ++i)
	{
		FESurface* surf = m_surf[i];
		double v = (double)i;

		for (int n = 0; n < surf->Nodes(); n++)
		{
			int nid = surf->NodeIndex(n);
			assert((nid >= 0) && (nid < NN));
			val_mesh[nid] = v;
			bn_mesh[nid] = 1;
		}
	}

	FENodeList nodeList = elset->GetNodeList();
	int nn = nodeList.Size();
	vector<int> bn(nn, 0);
	vector<double> val(nn, 0.0);
	for (int i = 0; i < nn; ++i)
	{
		int globalIndex = nodeList[i];
		bn[i] = bn_mesh[globalIndex];
		val[i] = val_mesh[globalIndex];
	}

	// count equations
	int neq = 0;
	for (int i = 0; i < nn; ++i)
	{
		if (bn[i] == 0)
		{
			bn[i] = neq++;
		}
		else
		{
			bn[i] = -neq-2;
			neq++;
		}
	}

	// solve Laplace equation
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	LinearSolver* ls = fecore.CreateDefaultLinearSolver(GetFEModel());
	if (ls == nullptr) return nullptr;

	SparseMatrix* A = ls->CreateSparseMatrix(Matrix_Type::REAL_SYMMETRIC);
	if (A == nullptr) return nullptr;

	FEGlobalMatrix K(A);
	FEModel dummy;
	std::vector<double> b(neq, 0);
	FELinearSystem LS(&dummy, K, b, val, true);

	// build matrix profile
	K.build_begin(neq);
	for (int i = 0; i < elset->Elements(); ++i)
	{
		FEElement& el = elset->Element(i);
		int ne = el.Nodes();
		std::vector<int> lm(ne, -1);
		for (int k=0; k<ne; ++k)
		{
			int l = nodeList.GlobalToLocalID(el.m_node[k]);
			lm[k] = bn[l];
		}
		K.build_add(lm);
	}
	K.build_end();
	A->Zero();

	// fill matrix and RHS
	FEElementMatrix ke;
	vector<int> lm;
	FEDomainList& domList = elset->GetDomainList();
	for (int ndom = 0; ndom < domList.size(); ++ndom)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(*domList.GetDomain(ndom));
		for (int m = 0; m < dom.Elements(); ++m)
		{
			// get the element
			FESolidElement& el = dom.Element(m);

			int neln = el.Nodes();

			// get the element stiffness matrix
			ke.resize(neln, neln);
			lm.resize(neln);
			vector<vec3d> gradN(neln);

			// calculate stiffness
			int nint = el.GaussPoints();

			// gauss weights
			double* w = el.GaussWeights();

			// nodal coordinates
			vec3d rt[FEElement::MAX_NODES];
			for (int j = 0; j < neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;

			// repeat over integration points
			ke.zero();
			for (int n = 0; n < nint; ++n)
			{
				double Jt = dom.ShapeGradient0(el, n, gradN.data());

				// calculate stiffness component
				for (int i = 0; i < neln; ++i) {
					for (int j = 0; j < neln; ++j)
						ke[i][j] += (gradN[i] * gradN[j]) * w[n] * Jt;
				}
			}

			// get the element's LM vector
			for (int j = 0; j < el.Nodes(); ++j)
			{
				int m = nodeList.GlobalToLocalID(el.m_node[j]);
				lm[j] = bn[m];
			}

			// assemble element matrix in global stiffness matrix
			ke.SetIndices(lm);
			LS.Assemble(ke);
		}
	}

	// solve linear system
	std::vector<double> x(neq, 0);
	ls->PreProcess();
	ls->Factor();
	ls->BackSolve(x, b);
	for (int i = 0; i < nn; ++i)
	{
		int m = bn[i];
		if (m >= 0) val[i] = x[m];
	}


	vec3d N = m_normal;
	N.unit();

	FEDomainMap* map = new FEDomainMap(FEDataType::FE_VEC3D, Storage_Fmt::FMT_ITEM);
	map->Create(elset);

	for (int i=0; i<elset->Elements(); ++i)
	{
		FESolidElement& el = dynamic_cast<FESolidElement&>(elset->Element(i));
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(*el.GetMeshPartition());

		int ne = el.Nodes();
		double vn[FEElement::MAX_NODES];
		for (int j = 0; j < ne; ++j)
		{
			int m = nodeList.GlobalToLocalID(el.m_node[j]);
			vn[j] = val[m];
		}

		// calculate average gradient
		vec3d grad(0,0,0);
		for (int j = 0; j < el.GaussPoints(); ++j)
		{
			grad += dom.gradient(el, vn, j);
		}
		grad.Normalize();

		if (m_cross)
		{
			vec3d b = grad ^ N;
			b.unit();
			grad = b;
		}

		// do plane rotations
		if (m_outAngle != 0)
		{
			vec3d a = grad ^ N;
			quatd q(m_outAngle * PI / 180.0, a);
			q.RotateVector(grad);
		}

		if (m_inAngle != 0)
		{
			quatd q(m_inAngle * PI / 180.0, N);
			q.RotateVector(grad);
		}

		map->setValue(i, grad);
	}
	return map;
}
