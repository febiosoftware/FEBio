#include "stdafx.h"
#include "Preconditioner.h"
#include "CompactUnSymmMatrix.h"

//-----------------------------------------------------------------------------
// We must undef PARDISO since it is defined as a function in mkl_solver.h
#ifdef MKL_ISS
#ifdef PARDISO
#undef PARDISO
#endif
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#endif // MKL_ISS

Preconditioner::Preconditioner()
{
}

Preconditioner::~Preconditioner()
{
}

//=================================================================================================

ILU0_Preconditioner::ILU0_Preconditioner()
{
	m_checkZeroDiagonal = true;
	m_zeroThreshold = 1e-16;
	m_zeroReplace = 1e-10;

	m_K = 0;
}

bool ILU0_Preconditioner::Create(SparseMatrix* A)
{
	m_K = dynamic_cast<CRSSparseMatrix*>(A);
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

void ILU0_Preconditioner::mult_vector(double* x, double* y)
{
	int ivar = m_K->Rows();
	int* ia = m_K->Pointers();
	int* ja = m_K->Indices();

	char cvar1 = 'L';
	char cvar = 'N';
	char cvar2 = 'U';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &m_bilu0[0], ia, ja, &x[0], &m_tmp[0]);
	cvar1 = 'U';
	cvar = 'N';
	cvar2 = 'N';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &m_bilu0[0], ia, ja, &m_tmp[0], &y[0]);
}

//=================================================================================================

ILUT_Preconditioner::ILUT_Preconditioner()
{
	m_maxfill = 1;
	m_fillTol = 1e-16;

	m_checkZeroDiagonal = true;
	m_zeroThreshold = 1e-16;
	m_zeroReplace = 1e-10;
}

bool ILUT_Preconditioner::Create(SparseMatrix* A)
{
	m_K = dynamic_cast<CRSSparseMatrix*>(A);
	if (m_K == 0) return false;

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
		dpar[31] = m_zeroReplace;
	}
	else
	{
		// do this to avoid a warning from the preconditioner
		dpar[30] = m_fillTol;
	}

	const int PCsize = (2 * m_maxfill + 1)*N + m_maxfill*(m_maxfill + 1) + 1;
	m_bilut.resize(PCsize, 0.0);
	m_jbilut.resize(PCsize, 0);
	m_ibilut.resize(N + 1, 0);

	m_tmp.resize(N, 0.0);

	int ierr;
	dcsrilut(&ivar, pa, ia, ja, &m_bilut[0], &m_ibilut[0], &m_jbilut[0], &m_fillTol, &m_maxfill, ipar, dpar, &ierr);
	if (ierr != 0) return false;

	return true;
}

void ILUT_Preconditioner::mult_vector(double* x, double* y)
{
	int ivar = m_K->Rows();
	char cvar1 = 'L';
	char cvar = 'N';
	char cvar2 = 'U';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &m_bilut[0], &m_ibilut[0], &m_jbilut[0], x, &m_tmp[0]);
	cvar1 = 'U';
	cvar = 'N';
	cvar2 = 'N';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &m_bilut[0], &m_ibilut[0], &m_jbilut[0], &m_tmp[0], y);
}


//=================================================================================================
DiagonalPreconditioner::DiagonalPreconditioner()
{
	m_P = 0;
}

// create a preconditioner for a sparse matrix
bool DiagonalPreconditioner::Create(SparseMatrix* A)
{
	if (A == 0) return false;

	int N = A->Rows();
	if (A->Columns() != N) return false;

	m_D.resize(N);
	for (int i=0; i<N; ++i)
	{
		double dii = A->diag(i);
		if (dii == 0.0) return false;
		m_D[i] = 1.0 / dii;
	}

	return true;
}

// apply to vector P x = y
void DiagonalPreconditioner::mult_vector(double* x, double* y)
{
	int N = (int)m_D.size();

#pragma omp parallel for
	for (int i=0; i<N; ++i)
	{
		y[i] = x[i]*m_D[i];
	}
}
