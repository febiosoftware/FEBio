// SparseMatrix.cpp: implementation of the SparseMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SparseMatrix.h"
#include <memory.h>

//-----------------------------------------------------------------------------
SparseMatrix::SparseMatrix()
{
	m_nrow = m_ncol = 0;
	m_nsize = 0;
}

SparseMatrix::~SparseMatrix()
{
}

void SparseMatrix::Clear()
{
	m_nrow = m_ncol = 0;
	m_nsize = 0;
}
