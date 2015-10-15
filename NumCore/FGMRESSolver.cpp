#include "stdafx.h"
#include "FGMRESSolver.h"

//-----------------------------------------------------------------------------
// We must undef PARDISO since it is defined as a function in mkl_solver.h
#ifdef MKL_ISS
#ifdef PARDISO
#undef PARDISO
#endif
#include "mkl_solver.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#endif // MKL_ISS

//-----------------------------------------------------------------------------
FGMRESSolver::FGMRESSolver() : m_pA(0)
{
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRESSolver::CreateSparseMatrix(Matrix_Type ntype)
{
#ifdef MKL_ISS
	if (ntype != REAL_UNSYMMETRIC) return 0;
	m_pA = new CompactUnSymmMatrix(1, true);
	return m_pA;
#else
	return 0;
#endif
}

//-----------------------------------------------------------------------------
bool FGMRESSolver::BackSolve(vector<double>& x, vector<double>& b)
{
#ifdef MKL_ISS
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// number of equations
	MKL_INT N = m_pA->Size();

	// data allocation
	MKL_INT ipar[128];
	double dpar[128];
	MKL_INT RCI_request;
	int M = (N < 150 ? N : 150); // this is the default value of par[15] (i.e. par[14] in C)
	vector<double> tmp(N*(2*M+1)+(M*(M+9))/2+1);
	double* ptmp = &tmp[0];
	MKL_INT ivar = N;

	// initialize the solver
	dfgmres_init(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// Set the desired parameters:
	ipar[8]=1;		// do residual stopping test
	ipar[9]=0;		// do not request for the user defined stopping test
	ipar[11]=1;		// do the check of the norm of the next generated vector automatically
//	dpar[0]=1.0E-3;	// set the relative tolerance to 1.0D-3 instead of default value 1.0D-6

	// Check the correctness and consistency of the newly set parameters
	dfgmres_check(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);
	if (RCI_request!=0) { MKL_Free_Buffers(); return false; }

	// solve the problem
	bool bdone = false;
	bool bconverged = false;
	while (!bdone)
	{
		// compute the solution via FGMRES
		dfgmres(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);

		switch (RCI_request)
		{
		case 0: // solution converged. 
			bdone = true;
			bconverged = true;
			break;
		case 1:
			{
				char cvar = 'N'; // multiply with unmodified A
				double* pa = m_pA->Values();
				int* ia = m_pA->Pointers();
				int* ja = m_pA->Indices();
				mkl_dcsrgemv(&cvar, &ivar, pa, ia, ja, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);
			}
			break;
		default:	// something went wrong
			bdone = true;
			bconverged = false;
		}
	}

	// get the solution. 
	MKL_INT itercount;
	dfgmres_get(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp, &itercount);

	MKL_Free_Buffers();
	return bconverged;

#else
	return false;
#endif // MKL_ISS
}

//=============================================================================
FGMRES_ILUT_Solver::FGMRES_ILUT_Solver() : m_pA(0)
{
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRES_ILUT_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
#ifdef MKL_ISS
	if (ntype != REAL_UNSYMMETRIC) return 0;
	m_pA = new CompactUnSymmMatrix(1, true);
	return m_pA;
#else
	return 0;
#endif
}

//-----------------------------------------------------------------------------
bool FGMRES_ILUT_Solver::BackSolve(vector<double>& x, vector<double>& b)
{
#ifdef MKL_ISS
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// number of equations
	MKL_INT N = m_pA->Size();
	double* pa = m_pA->Values();
	int* ia = m_pA->Pointers();
	int* ja = m_pA->Indices();

	// data allocation
	MKL_INT ipar[128];
	double dpar[128];
	MKL_INT RCI_request;
	int M = (N < 150 ? N : 150); // this is the default value of par[15] (i.e. par[14] in C)
	vector<double> tmp(N*(2*M+1)+(M*(M+9))/2+1);
	vector<double> trvec(N);
	double* ptmp = &tmp[0];
	MKL_INT ivar = N;

	// initialize the solver
	dfgmres_init(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// parameters affecting the pre-conditioner
	ipar[30]=1;			// change small diagonal value to that given by dpar[31]
	dpar[30]=1.E-5;		// override the default diagonal value that is used instead of zero

	// calculate the pre-conditioner
	int ierr;
	double tol = 1e-6;
	int maxfil = 1;
	const int PCsize = (2*maxfil+1)*N-maxfil*(maxfil+1)+1;
	vector<double> bilut(PCsize);
	vector<int> jbilut(PCsize);
	vector<int> ibilut(N+1);
	dcsrilut(&ivar, pa, ia, ja, &bilut[0], &ibilut[0], &jbilut[0], &tol, &maxfil, ipar, dpar, &ierr);
	if (ierr != 0) { MKL_Free_Buffers(); return false; }

	// Set the desired parameters:
	ipar[8]=1;		// do residual stopping test
	ipar[9]=0;		// do not request for the user defined stopping test
	ipar[10]=1;		// do the pre-conditioned version of the FGMRES iterative solver
	ipar[11]=1;		// do the check of the norm of the next generated vector automatically
//	dpar[0]=1.0E-3;	// set the relative tolerance to 1.0D-3 instead of default value 1.0D-6

	// Check the correctness and consistency of the newly set parameters
	dfgmres_check(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);
	if (RCI_request!=0) { MKL_Free_Buffers(); return false; }

	// solve the problem
	bool bdone = false;
	bool bconverged = false;
	while (!bdone)
	{
		// compute the solution via FGMRES
		dfgmres(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);

		switch (RCI_request)
		{
		case 0: // solution converged. 
			bdone = true;
			bconverged = true;
			break;
		case 1:
			{
				char cvar = 'N'; // multiply with unmodified A
				mkl_dcsrgemv(&cvar, &ivar, pa, ia, ja, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);
			}
			break;
		case 3:	// do the pre-conditioning step
			{
				char cvar1='L';
				char cvar='N';
				char cvar2='U';
				mkl_dcsrtrsv(&cvar1,&cvar,&cvar2,&ivar,&bilut[0],&ibilut[0],&jbilut[0],&tmp[ipar[21]-1], &trvec[0]);
				cvar1='U';
				cvar='N';
				cvar2='N';
				mkl_dcsrtrsv(&cvar1,&cvar,&cvar2,&ivar,&bilut[0],&ibilut[0],&jbilut[0], &trvec[0], &tmp[ipar[22]-1]);
			}
			break;
		default:	// something went wrong
			bdone = true;
			bconverged = false;
		}
	}

	// get the solution. 
	MKL_INT itercount;
	dfgmres_get(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp, &itercount);

	MKL_Free_Buffers();
	return bconverged;

#else
	return false;
#endif // MKL_ISS
}

//=============================================================================
FGMRES_ILU0_Solver::FGMRES_ILU0_Solver() : m_pA(0)
{
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRES_ILU0_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
#ifdef MKL_ISS
	if (ntype != REAL_UNSYMMETRIC) return 0;
	m_pA = new CompactUnSymmMatrix(1, true);
	return m_pA;
#else
	return 0;
#endif
}

//-----------------------------------------------------------------------------
bool FGMRES_ILU0_Solver::BackSolve(vector<double>& x, vector<double>& b)
{
#ifdef MKL_ISS
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// number of equations
	MKL_INT N = m_pA->Size();
	MKL_INT NNZ = m_pA->NonZeroes();
	double* pa = m_pA->Values();
	int* ia = m_pA->Pointers();
	int* ja = m_pA->Indices();

	// data allocation
	MKL_INT ipar[128];
	double dpar[128];
	MKL_INT RCI_request;
	int M = (N < 150 ? N : 150); // this is the default value of par[15] (i.e. par[14] in C)
	vector<double> tmp(N*(2*M+1)+(M*(M+9))/2+1);
	vector<double> trvec(N);
	double* ptmp = &tmp[0];
	MKL_INT ivar = N;

	// initialize the solver
	dfgmres_init(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// parameters affecting the pre-conditioner
	ipar[30]=1;			// change small diagonal value to that given by dpar[31]
	dpar[30]=1.E-20;		// override the default diagonal value that is used instead of zero
	dpar[31]=1.E-16;		// override the default diagonal value that is used instead of zero

	// calculate the pre-conditioner
	int ierr;
	vector<double> bilu0(NNZ);
	dcsrilu0(&ivar, pa, ia, ja, &bilu0[0], ipar, dpar, &ierr);
	if (ierr != 0) { MKL_Free_Buffers(); return false; }

	// Set the desired parameters:
	ipar[8]=1;		// do residual stopping test
	ipar[9]=0;		// do not request for the user defined stopping test
	ipar[10]=1;		// do the pre-conditioned version of the FGMRES iterative solver
	ipar[11]=1;		// do the check of the norm of the next generated vector automatically
//	dpar[0]=1.0E-3;	// set the relative tolerance to 1.0D-3 instead of default value 1.0D-6

	// Check the correctness and consistency of the newly set parameters
	dfgmres_check(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);
	if (RCI_request!=0) { MKL_Free_Buffers(); return false; }

	// solve the problem
	bool bdone = false;
	bool bconverged = false;
	while (!bdone)
	{
		// compute the solution via FGMRES
		dfgmres(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);

		switch (RCI_request)
		{
		case 0: // solution converged. 
			bdone = true;
			bconverged = true;
			break;
		case 1:
			{
				char cvar = 'N'; // multiply with unmodified A
				mkl_dcsrgemv(&cvar, &ivar, pa, ia, ja, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);
			}
			break;
		case 3:	// do the pre-conditioning step
			{
				char cvar1='L';
				char cvar='N';
				char cvar2='U';
				mkl_dcsrtrsv(&cvar1,&cvar,&cvar2,&ivar,&bilu0[0], ia, ja, &tmp[ipar[21]-1], &trvec[0]);
				cvar1='U';
				cvar='N';
				cvar2='N';
				mkl_dcsrtrsv(&cvar1,&cvar,&cvar2,&ivar,&bilu0[0], ia, ja, &trvec[0], &tmp[ipar[22]-1]);
			}
			break;
		default:	// something went wrong
			bdone = true;
			bconverged = false;
		}
	}

	// get the solution. 
	MKL_INT itercount;
	dfgmres_get(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp, &itercount);

	MKL_Free_Buffers();
	return bconverged;

#else
	return false;
#endif // MKL_ISS
}
