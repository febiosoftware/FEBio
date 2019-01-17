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

ILUT_Preconditioner::ILUT_Preconditioner(FEModel* fem) : Preconditioner(fem)
{
	m_maxfill = 1;
	m_fillTol = 1e-16;

	m_checkZeroDiagonal = true;
	m_zeroThreshold = 1e-16;
	m_zeroReplace = 1e-10;
}

bool ILUT_Preconditioner::Create()
{
	m_K = dynamic_cast<CRSSparseMatrix*>(GetSparseMatrix());
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

bool ILUT_Preconditioner::mult_vector(double* x, double* y)
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

	return true;
}
