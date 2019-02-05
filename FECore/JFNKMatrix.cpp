#include "stdafx.h"
#include "JFNKMatrix.h"
#include "FENewtonSolver.h"
#include "FEModel.h"

JFNKMatrix::JFNKMatrix(FENewtonSolver* pns, SparseMatrix* K) : m_pns(pns), m_K(K)
{
	m_nrow = m_ncol = pns->m_neq;
	m_nsize = 0;
	m_eps = pns->m_jfnk_eps;

	m_policy = ZERO_PRESCRIBED_DOFS;

	// TODO: For contact problems we'll need some mechanism to change the array size
	m_v.resize(m_nrow);
	m_R.resize(m_nrow);

	// figure out the free and prescribed equation numbers
	m_freeNodes.clear();
	m_prescribedNodes.clear();

	FEModel* fem = m_pns->GetFEModel();
	FEMesh& mesh = fem->GetMesh();
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j = 0; j < node.m_ID.size(); ++j)
		{
			int id = node.m_ID[j];

			if (id >= 0) m_freeNodes.push_back(id);
			if (id < -1) m_prescribedNodes.push_back(-id - 2);
		}
	}

	// make sure it all matches
	assert(m_freeNodes.size() + m_prescribedNodes.size() == m_pns->m_neq);
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
		for (int i = 0; i < m_freeNodes.size(); ++i)
		{
			int id = m_freeNodes[i];
			m_v[id] = m_eps*x[id];
		}
		for (int i = 0; i < m_prescribedNodes.size(); ++i)
		{
			int id = m_prescribedNodes[i];
			m_v[id] = 0.0;
		}
	}
	else
	{
		for (int i = 0; i < m_freeNodes.size(); ++i)
		{
			int id = m_freeNodes[i];
			m_v[id] = 0.0;
		}
		for (int i = 0; i < m_prescribedNodes.size(); ++i)
		{
			int id = m_prescribedNodes[i];
			m_v[id] = m_eps*x[id];
		}
	}

	m_pns->Update2(m_v);
	if (m_pns->Residual(m_R) == false) return false;

	for (int i = 0; i < m_freeNodes.size(); ++i)
	{
		int id = m_freeNodes[i];
		r[id] = (m_R0[id] - m_R[id]) / m_eps;
	}

	if (m_policy == ZERO_PRESCRIBED_DOFS)
	{
		for (int i = 0; i < m_prescribedNodes.size(); ++i)
		{
			int id = m_prescribedNodes[i];
			r[id] = x[id];
		}
	}
	else
	{
		for (int i = 0; i < m_prescribedNodes.size(); ++i)
		{
			int id = m_prescribedNodes[i];
			r[id] = 0.0;
		}
	}

	return  true;
}
