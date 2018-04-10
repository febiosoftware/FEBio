#include "stdafx.h"
#include "FENewtonStrategy.h"
#include "FENewtonSolver.h"

FENewtonStrategy::FENewtonStrategy(FENewtonSolver* pns) : m_pns(pns)
{
	m_maxups = 10;
	m_max_buf_size = 0; // when zero, it should default to m_maxups
	m_cycle_buffer = true;

	m_nups = 0;
}

FENewtonStrategy::~FENewtonStrategy()
{
}

//! initialize the linear system
SparseMatrix* FENewtonStrategy::CreateSparseMatrix(Matrix_Type mtype)
{
	if (m_pns == 0) return 0;

	LinearSolver* plinsolve = m_pns->m_plinsolve;

	SparseMatrix* pS = plinsolve->CreateSparseMatrix(mtype);

	return pS;
}

bool FENewtonStrategy::ReformStiffness()
{
	return m_pns->ReformStiffness();
}
