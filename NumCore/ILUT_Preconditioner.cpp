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
#include "ILUT_Preconditioner.h"
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

BEGIN_FECORE_CLASS(ILUT_Preconditioner, Preconditioner)
	ADD_PARAMETER(m_maxfill, "maxfill");
	ADD_PARAMETER(m_fillTol, "filltol");
	ADD_PARAMETER(m_checkZeroDiagonal, "replace_zero_diagonal");
	ADD_PARAMETER(m_zeroThreshold    , "zero_threshold");
	ADD_PARAMETER(m_zeroReplace      , "zero_replace");
END_FECORE_CLASS();

ILUT_Preconditioner::ILUT_Preconditioner(FEModel* fem) : Preconditioner(fem)
{
	m_maxfill = 1;
	m_fillTol = 1e-16;

	m_checkZeroDiagonal = true;
	m_zeroThreshold = 1e-16;
	m_zeroReplace = 1e-10;
}

SparseMatrix* ILUT_Preconditioner::CreateSparseMatrix(Matrix_Type ntype)
{
	if (ntype != REAL_UNSYMMETRIC) return nullptr;
	m_K = new CRSSparseMatrix(1);
	SetSparseMatrix(m_K);
	return m_K;
}

#ifdef MKL_ISS

bool ILUT_Preconditioner::Factor()
{
	if (m_K == 0) return false;
	assert(m_K->Offset() == 1);

	int N = m_K->Rows();
	int ivar = N;

	double* pa = m_K->Values();
	int* ia = m_K->Pointers();
	int* ja = m_K->Indices();

	MKL_INT ipar[128] = { 0 };
	double dpar[128] = { 0.0 };

	if (m_checkZeroDiagonal)
	{
		ipar[30] = 1;
		dpar[30] = m_zeroThreshold;
	}
	else
	{
		// do this to avoid a warning from the preconditioner
		dpar[30] = m_fillTol;
	}

	const int PCsize = (2 * m_maxfill + 1)*N - m_maxfill*(m_maxfill + 1) + 1;
	m_bilut.resize(PCsize, 0.0);
	m_jbilut.resize(PCsize, 0);
	m_ibilut.resize(N + 1, 0);

	m_tmp.resize(N, 0.0);

	int ierr;
	dcsrilut(&ivar, pa, ia, ja, &m_bilut[0], &m_ibilut[0], &m_jbilut[0], &m_fillTol, &m_maxfill, ipar, dpar, &ierr);
	if (ierr != 0) return false;

	return true;
}

bool ILUT_Preconditioner::BackSolve(double* x, double* y)
{
	int ivar = m_K->Rows();
	char cvar1 = 'L';
	char cvar = 'N';
	char cvar2 = 'U';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &m_bilut[0], &m_ibilut[0], &m_jbilut[0], y, &m_tmp[0]);
	cvar1 = 'U';
	cvar = 'N';
	cvar2 = 'N';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &m_bilut[0], &m_ibilut[0], &m_jbilut[0], &m_tmp[0], x);

	return true;
}
#else
bool ILUT_Preconditioner::Factor() { return false; }
bool ILUT_Preconditioner::BackSolve(double* x, double* y) { return false; }
#endif