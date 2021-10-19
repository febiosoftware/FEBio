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
#include "JFNKMatrix.h"
#include "FENewtonSolver.h"
#include "FEModel.h"
#include "FEDomain.h"

JFNKMatrix::JFNKMatrix(FENewtonSolver* pns, SparseMatrix* K) : m_pns(pns), m_K(K)
{
	m_nrow = m_ncol = pns->m_neq;
	m_nsize = 0;

	m_bauto_eps = false;
	m_eps = 1e-6;

	m_policy = ZERO_PRESCRIBED_DOFS;

	// TODO: For contact problems we'll need some mechanism to change the array size
	m_v.resize(m_nrow);
	m_R.resize(m_nrow);

	// figure out the free and prescribed equation numbers
	m_freeDofs.clear();
	m_prescribedDofs.clear();

	FEModel* fem = m_pns->GetFEModel();
	FEMesh& mesh = fem->GetMesh();
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j = 0; j < node.m_ID.size(); ++j)
		{
			int id = node.m_ID[j];

			if (id >= 0) m_freeDofs.push_back(id);
			if (id < -1) m_prescribedDofs.push_back(-id - 2);
		}
	}

	// Add element dofs
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NEL = dom.Elements();
		for (int j = 0; j < NEL; ++j)
		{
			FEElement& elj = dom.ElementRef(j);
			if (elj.m_lm >= 0) m_freeDofs.push_back(elj.m_lm);
		}
	}

	// make sure it all matches
	assert(m_freeDofs.size() + m_prescribedDofs.size() == m_pns->m_neq);
}

//! set matrix policy
void JFNKMatrix::SetPolicy(MultiplyPolicy p)
{
	m_policy = p;
}

//! Set the forward difference epsilon
void JFNKMatrix::SetEpsilon(double eps)
{
	m_eps = eps;
}

//! Create a sparse matrix from a sparse-matrix profile
void JFNKMatrix::Create(SparseMatrixProfile& MP) 
{ 
	m_K->Create(MP); 
	m_nrow = m_K->Rows();
	m_ncol = m_K->Columns();
	m_nsize = m_K->NonZeroes(); 
}

void JFNKMatrix::SetReferenceResidual(std::vector<double>& R0)
{
	m_R0 = R0;
}

bool JFNKMatrix::mult_vector(double* x, double* r)
{
	int neq = (int)m_pns->m_ui.size();

	if (m_policy == ZERO_PRESCRIBED_DOFS)
	{
		for (int i = 0; i < m_freeDofs.size(); ++i)
		{
			int id = m_freeDofs[i];
			m_v[id] = x[id];
		}
		for (int i = 0; i < m_prescribedDofs.size(); ++i)
		{
			int id = m_prescribedDofs[i];
			m_v[id] = 0.0;
		}
	}
	else
	{
		for (int i = 0; i < m_freeDofs.size(); ++i)
		{
			int id = m_freeDofs[i];
			m_v[id] = 0.0;
		}
		for (int i = 0; i < m_prescribedDofs.size(); ++i)
		{
			int id = m_prescribedDofs[i];
			m_v[id] = x[id];
		}
	}

	double eps = m_eps;
	if (m_bauto_eps)
	{
		vector<double> u;
		m_pns->GetSolutionVector(u);

		assert(neq == Rows());
		double norm_v = l2_norm(m_v);

		eps = 0.0;
		if (norm_v != 0.0)
		{
			for (int i = 0; i < neq; ++i) eps += fabs(u[i]);
			eps *= m_eps / (neq*norm_v);
		}

		eps += m_eps;
	}

	// multiply by eps
	for (int i = 0; i < neq; ++i) m_v[i] *= eps;

	m_pns->Update2(m_v);
	if (m_pns->Residual(m_R) == false) return false;

	for (int i = 0; i < m_freeDofs.size(); ++i)
	{
		int id = m_freeDofs[i];
		r[id] = (m_R0[id] - m_R[id]) / eps;
	}

	if (m_policy == ZERO_PRESCRIBED_DOFS)
	{
		for (int i = 0; i < m_prescribedDofs.size(); ++i)
		{
			int id = m_prescribedDofs[i];
			r[id] = x[id];
		}
	}
	else
	{
		for (int i = 0; i < m_prescribedDofs.size(); ++i)
		{
			int id = m_prescribedDofs[i];
			r[id] = 0.0;
		}
	}

	return  true;
}
