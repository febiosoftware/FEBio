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
#include "ILU0_Preconditioner.h"
#include <FECore/CompactUnSymmMatrix.h>

// We must undef PARDISO since it is defined as a function in mkl_solver.h
#ifdef MKL_ISS
#ifdef PARDISO
#undef PARDISO
#endif
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#endif // MKL_ISS

BEGIN_FECORE_CLASS(ILU0_Preconditioner, Preconditioner)
	ADD_PARAMETER(m_checkZeroDiagonal, "replace_zero_diagonal");
	ADD_PARAMETER(m_zeroThreshold    , "zero_threshold");
	ADD_PARAMETER(m_zeroReplace      , "zero_replace");
END_FECORE_CLASS();

//=================================================================================================

ILU0_Preconditioner::ILU0_Preconditioner(FEModel* fem) : Preconditioner(fem)
{
	m_checkZeroDiagonal = true;
	m_zeroThreshold = 1e-16;
	m_zeroReplace = 1e-10;

	m_K = 0;
}

SparseMatrix* ILU0_Preconditioner::CreateSparseMatrix(Matrix_Type ntype)
{
	if (ntype != REAL_UNSYMMETRIC) return nullptr;
	m_K = new CRSSparseMatrix(1);
	return m_K;
}

#ifdef MKL_ISS
bool ILU0_Preconditioner::Factor()
{

	if (m_K == 0) return false;
	assert(m_K->Offset() == 1);

	int N = m_K->Rows();
	int NNZ = m_K->NonZeroes();

	double* pa = m_K->Values();
	int* ia = m_K->Pointers();
	int* ja = m_K->Indices();

	MKL_INT ipar[128] = { 0 };
	double dpar[128] = { 0.0 };

	// parameters affecting the pre-conditioner
	if (m_checkZeroDiagonal)
	{
		ipar[30] = 1;
		dpar[30] = m_zeroThreshold;
		dpar[31] = m_zeroReplace;
	}

	m_tmp.resize(N, 0.0);

	m_bilu0.resize(NNZ);
	int ierr = 0;
	dcsrilu0(&N, pa, ia, ja, &m_bilu0[0], ipar, dpar, &ierr);
	if (ierr != 0) return false;

	return true;
} 

bool ILU0_Preconditioner::BackSolve(double* x, double* y)
{
	int ivar = m_K->Rows();
	int* ia = m_K->Pointers();
	int* ja = m_K->Indices();

	char cvar1 = 'L';
	char cvar = 'N';
	char cvar2 = 'U';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &m_bilu0[0], ia, ja, &y[0], &m_tmp[0]);
	cvar1 = 'U';
	cvar = 'N';
	cvar2 = 'N';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &m_bilu0[0], ia, ja, &m_tmp[0], &x[0]);

	return true;
}

#else
bool ILU0_Preconditioner::Factor() { return false; }
bool ILU0_Preconditioner::BackSolve(double* x, double* y) { return false; }
#endif
