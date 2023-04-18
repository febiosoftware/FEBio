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
#include "FGMRESSolver.h"
#include <FECore/CompactSymmMatrix.h>
#include <FECore/CompactUnSymmMatrix.h>
#include <FECore/log.h>
#include "MatrixTools.h"

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

BEGIN_FECORE_CLASS(FGMRESSolver, IterativeLinearSolver)
	ADD_PARAMETER(m_maxiter       , "max_iter");
	ADD_PARAMETER(m_print_level   , "print_level");
	ADD_PARAMETER(m_doResidualTest, "check_residual");
	ADD_PARAMETER(m_nrestart      , "max_restart");
	ADD_PARAMETER(m_reltol        , "tol");
	ADD_PARAMETER(m_abstol        , "abs_tol");
	ADD_PARAMETER(m_maxIterFail   , "fail_max_iters");

	ADD_PROPERTY(m_P, "pc_left")->SetFlags(FEProperty::Optional);
	ADD_PROPERTY(m_R, "pc_right")->SetFlags(FEProperty::Optional);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FGMRESSolver::FGMRESSolver(FEModel* fem) : IterativeLinearSolver(fem), m_pA(0)
{
	m_maxiter = 0; // use default min(N, 150)
	m_print_level = 0;
	m_doResidualTest = true;
	m_doZeroNormTest = true;
	m_reltol = 0.0;
	m_abstol = 0.0;
	m_nrestart = 0; // use default = maxiter
	m_print_cn = false;

	m_do_jacobi = false;

	m_P = 0;	// we don't use a preconditioner for this solver
	m_R = 0;	// no right preconditioner

	m_maxIterFail = true;
}

//-----------------------------------------------------------------------------
// set the preconditioner
void FGMRESSolver::SetLeftPreconditioner(LinearSolver* P)
{
	m_P = P;
}

//-----------------------------------------------------------------------------
//! Set the right preconditioner
void FGMRESSolver::SetRightPreconditioner(LinearSolver* R)
{
	m_R = R;
}

//-----------------------------------------------------------------------------
// get the preconditioner
LinearSolver* FGMRESSolver::GetLeftPreconditioner()
{
	return m_P;
}

//-----------------------------------------------------------------------------
// get the preconditioner
LinearSolver* FGMRESSolver::GetRightPreconditioner()
{
	return m_R;
}

//-----------------------------------------------------------------------------
//! Set max nr of iterations
void FGMRESSolver::SetMaxIterations(int n)
{
	m_maxiter = n;
}

//-----------------------------------------------------------------------------
//! Get the max nr of iterations
int FGMRESSolver::GetMaxIterations() const
{
	return m_maxiter;
}

//-----------------------------------------------------------------------------
//! Set the nr of non-restarted iterations
void FGMRESSolver::SetNonRestartedIterations(int n)
{
	m_nrestart = n;
}

//-----------------------------------------------------------------------------
// Set the print level
void FGMRESSolver::SetPrintLevel(int n)
{
	m_print_level = n;
}

//-----------------------------------------------------------------------------
// set residual stopping test flag
void FGMRESSolver::DoResidualStoppingTest(bool b)
{
	m_doResidualTest = b;
}

//-----------------------------------------------------------------------------
// set zero norm stopping test flag
void FGMRESSolver::DoZeroNormStoppingTest(bool b)
{
	m_doZeroNormTest = b;
}

//-----------------------------------------------------------------------------
// set the convergence tolerance for the residual stopping test
void FGMRESSolver::SetRelativeResidualTolerance(double tol)
{
	m_reltol = tol;
}

//-----------------------------------------------------------------------------
// set the absolute convergence tolerance for the residual stopping test
void FGMRESSolver::SetAbsoluteResidualTolerance(double tol)
{
	m_abstol = tol;
}

//-----------------------------------------------------------------------------
//! This solver does not use a preconditioner
bool FGMRESSolver::HasPreconditioner() const 
{
	return ((m_P != 0) || (m_R != 0)); 
}

//-----------------------------------------------------------------------------
void FGMRESSolver::FailOnMaxIterations(bool b)
{
	m_maxIterFail = b;
}

//-----------------------------------------------------------------------------
void FGMRESSolver::PrintConditionNumber(bool b)
{
	m_print_cn = b;
}

//-----------------------------------------------------------------------------
// do jacobi preconditioning
void FGMRESSolver::DoJacobiPreconditioning(bool b)
{
	m_do_jacobi = b;
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRESSolver::CreateSparseMatrix(Matrix_Type ntype)
{
#ifdef MKL_ISS
	// Cleanup if necessary
	if (m_pA) delete m_pA; 
	m_pA = nullptr;

	// since FMGRES doesn't really care what matrix is requested, 
	// see if the preconditioner cares.
	if (m_P)
	{
		m_P->SetPartitions(m_part);
		m_pA = m_P->CreateSparseMatrix(ntype);
		return m_pA;
	}
	else if (m_R)
	{
		m_R->SetPartitions(m_part);
		m_pA = m_R->CreateSparseMatrix(ntype);
		return m_pA;
	}

	// if the matrix is still zero, let's just allocate one
	if (m_pA == nullptr)
	{
		// allocate new matrix
		switch (ntype)
		{
		case REAL_SYMMETRIC: m_pA = new CompactSymmMatrix(1); break;
		case REAL_UNSYMMETRIC: m_pA = new CRSSparseMatrix(1); break;
		case REAL_SYMM_STRUCTURE: m_pA = new CRSSparseMatrix(1); break;
		}
	}

	// return the matrix (Can be null if matrix format not supported!)
	return m_pA;
#else
	return 0;
#endif
}

//-----------------------------------------------------------------------------
bool FGMRESSolver::SetSparseMatrix(SparseMatrix* pA)
{
	m_pA = pA;
	return (m_pA != 0);
}

//-----------------------------------------------------------------------------
//! Clean up
void FGMRESSolver::Destroy()
{
	m_tmp.clear();
	m_tmp.shrink_to_fit();
}

//-----------------------------------------------------------------------------
bool FGMRESSolver::PreProcess() 
{
#ifdef MKL_ISS
	// number of equations
	MKL_INT N = m_pA->Rows();

	int M = (N < 150 ? N : 150); // this is the default value of ipar[14]

	if (m_nrestart > 0) M = m_nrestart;
	else if (m_maxiter > 0) M = m_maxiter;

	// allocate temp storage
	m_tmp.resize((N*(2 * M + 1) + (M*(M + 9)) / 2 + 1));

	m_Rv.resize(N);

	m_W.resize(N, 1.0);

	return true; 
#else
	return false;
#endif
}


//! Factor the matrix
bool FGMRESSolver::Factor()
{
	int neq = m_pA->Rows();
	if (m_do_jacobi)
	{
		for (int i = 0; i < neq; ++i)
		{
			double dii = fabs(m_pA->diag(i));
			if (dii == 0.0) return false;
			m_W[i] = 1.0 / sqrt(dii);
		}
	}

	if (m_do_jacobi)
		m_pA->scale(m_W, m_W);

	if (m_print_cn)
	{
		double c = NumCore::estimateConditionNumber(GetSparseMatrix());
		feLog("\tcondition number (est.) ................... : %lg\n\n", c);
	}

	// call the preconditioner
	if (m_P)
	{
		m_P->SetFEModel(GetFEModel());
		if (m_P->PreProcess() == false) return false;
		if (m_P->Factor() == false) return false;
	}

	if (m_R)
	{
		m_R->SetFEModel(GetFEModel());
		if (m_R->PreProcess() == false) return false;
		if (m_R->Factor() == false) return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FGMRESSolver::BackSolve(double* x, double* b)
{
#ifdef MKL_ISS
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// number of equations
	MKL_INT N = m_pA->Rows();

	// data allocation
	int M = (N < 150 ? N : 150); // this is the default value of ipar[4] and ipar[14]

	int nrestart = M;
	if (m_nrestart > 0) nrestart = m_nrestart;
	else if (m_maxiter > 0) nrestart = m_maxiter;

	int maxIter = M;
	if (m_maxiter > 0) maxIter = m_maxiter;

	// scale rhs
	vector<double> F(N);
	for (int i = 0; i < N; ++i) F[i] = m_W[i] * b[i];

	// initialize the solver
	MKL_INT ipar[128] = { 0 };
	double dpar[128] = { 0.0 };
	MKL_INT ivar = N;
	MKL_INT RCI_request;
	dfgmres_init(&ivar, &x[0], &F[0], &RCI_request, ipar, dpar, &m_tmp[0]);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// Set the desired parameters:
	ipar[ 4] = maxIter;	                        // max number of iterations
	ipar[ 7] = 1;								// do the stopping test for maximal number of iterations
	ipar[ 8] = (m_doResidualTest ? 1 : 0);		// do residual stopping test
	ipar[ 9] = 0;								// do not request for the user defined stopping test
	ipar[10] = (m_P != 0 ? 1 : 0);				// do the pre-conditioned version of the FGMRES iterative solver
	ipar[11] = (m_doZeroNormTest ? 1 : 0);		// do the check of the norm of the next generated vector automatically
	ipar[14] = nrestart;	                    // number of non-restarted iterations
	if (m_reltol > 0) dpar[0] = m_reltol;		// set the relative tolerance
	if (m_abstol > 0) dpar[1] = m_abstol;		// set the absolute tolerance

	// Check the correctness and consistency of the newly set parameters
	dfgmres_check(&ivar, &x[0], &F[0], &RCI_request, ipar, dpar, &m_tmp[0]);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// zero solution vector
	for (int i = 0; i < N; ++i) x[i] = 0.0;

	if (m_print_level > 0) feLog("FGMRES:\n");

	// solve the problem
	bool bdone = false;
	bool bconverged = !m_maxIterFail;
	while (!bdone)
	{
		// compute the solution via FGMRES
		dfgmres(&ivar, &x[0], &F[0], &RCI_request, ipar, dpar, &m_tmp[0]);

		switch (RCI_request)
		{
		case 0: // solution converged. 
			bdone = true;
			bconverged = true;
			break;
		case 1:
			{
				// do matrix-vector multiplication
				if (m_R)
				{
					// first apply the right preconditioner
					m_R->mult_vector(&m_tmp[ipar[21] - 1], &m_Rv[0]);

					// then multiply with matrix
					m_pA->mult_vector(&m_Rv[0], &m_tmp[ipar[22] - 1]);
				}
				else m_pA->mult_vector(&m_tmp[ipar[21] - 1], &m_tmp[ipar[22] - 1]);

				if (m_print_level > 1)
				{
					feLog("%3d = %lg (%lg), %lg (%lg)\n", ipar[3], dpar[4], dpar[3], dpar[6], dpar[7]);
				}
			}
			break;
		case 3:	// do the pre-conditioning step
			{
				assert(m_P);
				if (m_P->mult_vector(&m_tmp[ipar[21] - 1], &m_tmp[ipar[22] - 1]) == false)
				{
					bdone = true;
					bconverged = false;
				}
			}
			break;
		case 4:
			break;
		default:	// something went wrong
			bdone = true;
			bconverged = !m_maxIterFail;
		}
	}

	// get the solution. 
	MKL_INT itercount;

	dfgmres_get(&ivar, &x[0], &F[0], &RCI_request, ipar, dpar, &m_tmp[0], &itercount);

	if (m_do_jacobi)
	{
		for (int i = 0; i < N; ++i) x[i] *= m_W[i];
	}

	if (m_R)
	{
		m_R->mult_vector(&x[0], &m_Rv[0]);
		for (int i = 0; i < N; ++i) x[i] = m_Rv[i];
	}

	if (m_print_level > 0)
	{
		feLog("%3d = %lg (%lg), %lg (%lg)\n", ipar[3]+1, dpar[4], dpar[3], dpar[6], dpar[7]);
	}

//	MKL_Free_Buffers();

	// update stats
	UpdateStats(itercount);

	return bconverged;

#else
	return false;
#endif // MKL_ISS
}

//! convenience function for solving linear system Ax = b
bool FGMRESSolver::Solve(SparseMatrix* A, vector<double>& x, vector<double>& b)
{
	SetSparseMatrix(A);
	if (PreProcess() == false) return false;
	if (Factor() == false) return false;
	return BackSolve(&x[0], &b[0]);
}
