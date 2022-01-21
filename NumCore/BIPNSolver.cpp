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
#include "BIPNSolver.h"
#include <FECore/vector.h>
#include "RCICGSolver.h"
#include <FECore/SchurComplement.h>
#include <FECore/log.h>
#include "ILU0_Preconditioner.h"
#include "IncompleteCholesky.h"
#include "FGMRESSolver.h"
#include "BoomerAMGSolver.h"
#include <FECore/Preconditioner.h>

#ifdef MKL_ISS

// We must undef PARDISO since it is defined as a function in mkl_solver.h
#ifdef PARDISO
#undef PARDISO
#endif

#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"

// in SchurSolver.cpp
bool BuildDiagonalMassMatrix(FEModel* fem, BlockMatrix* K, CompactSymmMatrix* M, double scale);
bool BuildMassMatrix(FEModel* fem, BlockMatrix* K, CompactSymmMatrix* M, double scale);

BEGIN_FECORE_CLASS(BIPNSolver, LinearSolver)
	ADD_PARAMETER(m_maxiter             , "maxiter"    );
	ADD_PARAMETER(m_tol                 , "tol"        );
	ADD_PARAMETER(m_print_level         , "print_level");
	ADD_PARAMETER(m_use_cg              , "use_cg");
	ADD_PARAMETER(m_cg_maxiter          , "cg_maxiter" );
	ADD_PARAMETER(m_cg_tol              , "cg_tol"     );
	ADD_PARAMETER(m_cg_doResidualTest   , "cg_check_residual");
	ADD_PARAMETER(m_gmres_maxiter       , "gmres_maxiter");
	ADD_PARAMETER(m_gmres_tol           , "gmres_tol" );
	ADD_PARAMETER(m_gmres_doResidualTest, "gmres_check_residual");
	ADD_PARAMETER(m_gmres_pc            , "gmres_precondition");
	ADD_PARAMETER(m_do_jacobi           , "do_jacobi");
	ADD_PARAMETER(m_precondition_schur  , "precondition_schur");
END_FECORE_CLASS();

// constructor
BIPNSolver::BIPNSolver(FEModel* fem) : LinearSolver(fem), m_A(0)
{
	m_print_level = 0;
	m_maxiter = 10;
	m_tol = 1e-6;

	m_cg_maxiter = 0;
	m_cg_tol = 0.0;
	m_cg_doResidualTest = true;

	m_gmres_maxiter = 0;
	m_gmres_tol = 0.0;
	m_gmres_doResidualTest = true;
	m_gmres_pc = 0;

	m_do_jacobi = true;

	m_precondition_schur = 0;

	m_use_cg = true;

	m_Asolver = nullptr;
	m_PS = nullptr;
}

// set the max nr of BIPN iterations
void BIPNSolver::SetMaxIterations(int n)
{
	m_maxiter = n;
}

// set the output level
void BIPNSolver::SetPrintLevel(int n)
{
	m_print_level = n;
}

// Set the BIPN convergence tolerance
void BIPNSolver::SetTolerance(double eps)
{
	m_tol = eps;
}

// Use CG for step 2 or not
void BIPNSolver::UseConjugateGradient(bool b)
{
	m_use_cg = b;
}

// set the CG convergence parameters
void BIPNSolver::SetCGParameters(int maxiter, double tolerance, bool doResidualStoppingTest)
{
	m_cg_maxiter        = maxiter;
	m_cg_tol            = tolerance;
	m_cg_doResidualTest = doResidualStoppingTest;
}

// set the GMRES convergence parameters
void BIPNSolver::SetGMRESParameters(int maxiter, double tolerance, bool doResidualStoppingTest, int precondition)
{
	m_gmres_maxiter        = maxiter;
	m_gmres_tol            = tolerance;
	m_gmres_doResidualTest = doResidualStoppingTest;
	m_gmres_pc             = precondition;
}

// Do Jacobi preconditioner
void BIPNSolver::DoJacobiPreconditioner(bool b)
{
	m_do_jacobi = b;
}

// set the schur preconditioner option
void BIPNSolver::SetSchurPreconditioner(int n)
{
	m_precondition_schur = n;
}

//! Return a sparse matrix compatible with this solver
SparseMatrix* BIPNSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	// make sure we can support this matrix
	if (ntype != Matrix_Type::REAL_UNSYMMETRIC) return 0;

	// make sure we have two partitions
	if (m_part.size() != 2) return 0;

	int noffset = 1;
	if (m_gmres_pc == 2) noffset = 0;

	// allocate new matrix
	if (m_A) delete m_A;
	m_A = new BlockMatrix();
	m_A->Partition(m_part, ntype, noffset);

	// and return
	return m_A;
}

// set the sparse matrix
bool BIPNSolver::SetSparseMatrix(SparseMatrix* A)
{
	m_A = dynamic_cast<BlockMatrix*>(A);
	if (m_A == nullptr) return false;

	return true;
}

// allocate storage
bool BIPNSolver::PreProcess()
{
	// make sure we have a matrix
	if (m_A == 0) return false;
	BlockMatrix& A = *m_A;

	// get the number of equations
	int N = A.Rows();

	// make sure the partition is valid
	if (m_part.size() != 2) return false;
	int Nu = m_part[0];
	int Np = m_part[1];
	assert((Nu + Np) == N);

	// allocate pre-conditioners
	Kd.resize(Nu);
	Wm.resize(Nu);
	Wc.resize(Np);

	// allocate temp storage
	yu.resize(Nu);
	yp.resize(Np);

	yu_n.resize(Nu);
	yp_n.resize(Np);

	Rm0.resize(Nu);
	Rc0.resize(Np);

	Rm_n.resize(Nu);
	Rc_n.resize(Np);

	Yu.resize(m_maxiter, std::vector<double>(Nu));
	Yp.resize(m_maxiter, std::vector<double>(Np));

	RM.resize(Nu);
	RC.resize(Np);

	Rmu.resize(m_maxiter, std::vector<double>(Nu));
	Rmp.resize(m_maxiter, std::vector<double>(Nu));
	Rcu.resize(m_maxiter, std::vector<double>(Np));
	Rcp.resize(m_maxiter, std::vector<double>(Np));

	au.resize(m_maxiter+1);
	ap.resize(m_maxiter+1);

	du.resize(Nu);
	dp.resize(Np);

	// CG temp buffers
	if (m_use_cg)
		cg_tmp.resize(4* Np);
	else
	{
		// GMRES buffer
		int M = (Np < 150 ? Np : 150); // this is the default value of par[15] (i.e. par[14] in C)
		if (m_gmres_maxiter > 0) M = m_gmres_maxiter;

		// allocate temp storage
		cg_tmp.resize((Np*(2 * M + 1) + (M*(M + 9)) / 2 + 1));
	}

	// GMRES buffer
	int M = (Nu < 150 ? Nu : 150); // this is the default value of par[15] (i.e. par[14] in C)
	if (m_gmres_maxiter > 0) M = m_gmres_maxiter;

	// allocate temp storage
	gmres_tmp.resize((Nu*(2 * M + 1) + (M*(M + 9)) / 2 + 1));

	m_Asolver = new FGMRESSolver(GetFEModel());

	// initialize solver for A block
	LinearSolver* PC = nullptr;
	switch (m_gmres_pc)
	{
	case 0: PC = nullptr; break;
	case 1: PC = new ILU0_Preconditioner(GetFEModel()); break;
	case 2: PC = new BoomerAMGSolver(GetFEModel()); break;
	default:
		return false;
	}
		
	m_Asolver->SetMaxIterations(m_gmres_maxiter);
	m_Asolver->SetRelativeResidualTolerance(m_gmres_tol);
	if (m_Asolver->SetSparseMatrix(A.Block(0, 0).pA) == false) return false;

	if (m_Asolver->GetLeftPreconditioner())
	{
		m_Asolver->GetLeftPreconditioner()->SetSparseMatrix(A.Block(0, 0).pA);
	}
	if (m_Asolver->PreProcess() == false) return false;

	// preconditioner of Schur solver
	if ((m_precondition_schur > 0) && (m_PS == nullptr))
	{
		if (m_precondition_schur == 1)
		{
			// diagonal mass matrix
			CompactSymmMatrix* M = new CompactSymmMatrix(1);
			if (BuildDiagonalMassMatrix(GetFEModel(), m_A, M, 1.0) == false) return false;

			DiagonalPreconditioner* PS = new DiagonalPreconditioner(GetFEModel());
			if (PS->SetSparseMatrix(M) == false) return false;
			if (PS->PreProcess() == false) return false;
			if (PS->Factor() == false) return false;

			m_PS = PS;
		}
		else
		{
			// mass matrix
			CompactSymmMatrix* M = new CompactSymmMatrix(1);
			if (BuildMassMatrix(GetFEModel(), m_A, M, 1.0) == false) return false;

			// We do an incomplete cholesky factorization
			IncompleteCholesky* PS = new IncompleteCholesky(GetFEModel());
			if (PS->SetSparseMatrix(M) == false) return false;
			if (PS->PreProcess() == false) return false;
			if (PS->Factor() == false) return false;

			m_PS = PS;
		}
	}

	return true;
}

//! Pre-condition the matrix
bool BIPNSolver::Factor()
{
	// make sure we have a matrix
	if (m_A == 0) return false;
	BlockMatrix& A = *m_A;

	// get the number of equations
	int N = A.Rows();
	int Nu = m_part[0];
	int Np = m_part[1];

	// get the blocks
	// 
	//     | K | G |
	// A = |---+---|
	//     | D | L |
	//
	BlockMatrix::BLOCK& K = m_A->Block(0, 0);
	BlockMatrix::BLOCK& G = m_A->Block(0, 1);
	BlockMatrix::BLOCK& D = m_A->Block(1, 0);
	BlockMatrix::BLOCK& L = m_A->Block(1, 1);

	// calculate preconditioner
	Wm.assign(Nu, 1.0);
	Wc.assign(Np, 1.0);
	if (m_do_jacobi)
	{
		for (int i = 0; i < Nu; ++i)
		{
			double ki = fabs(K.pA->diag(i));
			if (ki > 0.0) 
				Wm[i] = 1.0 / sqrt(ki);
		}
		if (L.pA)
		{
			for (int i = 0; i < Np; ++i)
			{
				double li = fabs(L.pA->diag(i));
				if (li > 0.0)
					Wc[i] = 1.0 / sqrt(li);
			}
		}
	}

	// normalize the matrix
	K.pA->scale(Wm, Wm);
	G.pA->scale(Wm, Wc);
	D.pA->scale(Wc, Wm);
	if (L.pA) L.pA->scale(Wc, Wc);

	// Now calculate diagonal matrix of K
	if (m_do_jacobi) Kd.assign(Nu, 1.0);
	else
	{
		for (int i = 0; i < Nu; ++i)
		{
			double ki = K.pA->diag(i);
			if (ki == 0.0) return false;
			Kd[i] = 1.0 / ki;
		}
	}

	if (m_Asolver->Factor() == false) return false;

	return true;
}

//! Calculate the solution of RHS b and store solution in x
bool BIPNSolver::BackSolve(double* x, double* b)
{
	// make sure we have a matrix
	if (m_A == 0) return false;
	BlockMatrix& A = *m_A;

	BlockMatrix::BLOCK& K = m_A->Block(0, 0);
	BlockMatrix::BLOCK& G = m_A->Block(0, 1);
	BlockMatrix::BLOCK& D = m_A->Block(1, 0);
	BlockMatrix::BLOCK& L = m_A->Block(1, 1);

	// number of equations
	int N = A.Rows();
	int Nu = m_part[0];
	int Np = m_part[1];

	// normalize RHS
	for (int i = 0; i<Nu; ++i) Rm0[i] = Wm[i]*b[i     ];
	for (int i = 0; i<Np; ++i) Rc0[i] = Wc[i]*b[i + Nu];

	// initialize
	RM = Rm0;
	RC = Rc0;

	// calculate initial error
	double err_0 = Rm0*Rm0 + Rc0*Rc0, err_n = 0.0;

	// setup the Schur complement
	DiagonalPreconditioner DPC(nullptr);
	if (DPC.Create(K.pA) == false) return false;
	SchurComplementA S(&DPC, G.pA, D.pA, L.pA);
	S.NegateSchur(true);

	if (m_print_level != 0) feLog("--- Starting BIPN:\n");

	int MI = m_maxiter;
	int M = 2 * m_maxiter;
	matrix QM(M, M);
	QM.zero();

	// do the BIPN iterations 
	int niter = 0;
	for (int n=0; n<m_maxiter; ++n)
	{
		niter++;
		if (m_print_level != 0) feLog("BIPN %d: ", niter);

		// solve for yu_n (use GMRES): K*yu_n = RM[n]
		m_Asolver->ResetStats();
		m_Asolver->BackSolve(&yu_n[0], &(RM[0]));
		m_gmres1_iters = m_Asolver->GetStats().iterations;
		if (m_print_level != 0) feLog("%d, ", m_gmres1_iters);

		// compute the corrected residual Rc_n = RC[n] - D*yu_n
		D.vmult(yu_n, Rc_n);
		vsub(Rc_n, RC, Rc_n);

		// solve for yp_n (use CG): (L - D*G) * yp_n = Rc_n
		if (m_use_cg)
			m_cg_iters = cgsolve(&S, m_PS, yp_n, Rc_n);
		else
			m_cg_iters = gmressolve(&S, m_PS, yp_n, Rc_n);
		if (m_print_level != 0) feLog("%d, ", m_cg_iters);

		// compute corrected residual: Rm_n = RM[n] - G*yp_n
		G.vmult(yp_n, Rm_n);
		vsub(Rm_n, RM, Rm_n);

		// solve for yu_n (use GMRES): K*yu_n = Rm_n
		m_Asolver->ResetStats();
		m_Asolver->BackSolve(&yu_n[0], &Rm_n[0]);
		m_gmres2_iters = m_Asolver->GetStats().iterations;
		if (m_print_level != 0) feLog("%d, ", m_gmres2_iters);

		// calculate temp vectors
		K.vmult(yu_n, Rmu[n]);	// Rmu[n] = K*yu_n;
		G.vmult(yp_n, Rmp[n]);	// Rmp[n] = G*yp_n;
		D.vmult(yu_n, Rcu[n]);	// Rcu[n] = D*yu_n;
		L.vmult(yp_n, Rcp[n]);	// Rcp[n] = L*yp_n;

		// store solution candidates
		Yu[n] = yu_n;
		Yp[n] = yp_n;

		// update QM
		for (int j = 0; j <= n; ++j)
		{
			QM(n     , j     ) = Rmu[n] * Rmu[j] + Rcu[n] * Rcu[j];
			QM(n     , j + MI) = Rmu[n] * Rmp[j] + Rcu[n] * Rcp[j];
			QM(n + MI, j     ) = Rmp[n] * Rmu[j] + Rcp[n] * Rcu[j];
			QM(n + MI, j + MI) = Rmp[n] * Rmp[j] + Rcp[n] * Rcp[j];

			QM(j     , n     ) = QM(     n, j     );
			QM(j + MI, n     ) = QM(     n, j + MI);
			QM(j     , n + MI) = QM(n + MI, j     );
			QM(j + MI, n + MI) = QM(n + MI, j + MI);
		}

		// number of coefficients to determine
		const int m = 2 * (n + 1);

		// setup Q matrix
		matrix Q(m, m);
		vector<double> q(m);

		// fill in the blocks of Q and q
		for (int i = 0; i <= n; ++i)
		{
			for (int j = 0; j <= n; ++j)
			{
				Q(i        , j        ) = QM(i, j);
				Q(i        , j + n + 1) = QM(i, j+ MI);
				Q(i + n + 1, j        ) = QM(i+ MI, j);
				Q(i + n + 1, j + n + 1) = QM(i+ MI, j+ MI);
			}

			q[i        ] = Rm0*Rmu[i] + Rc0*Rcu[i];
			q[i + n + 1] = Rm0*Rmp[i] + Rc0*Rcp[i];
		}

		// solve for the coefficients
		vector<double> a(m);
		Q.solve(a, q);
		for (int i = 0; i <= n; ++i)
		{
			au[i] = a[i];
			ap[i] = a[i + n + 1];
		}

		// calculate error
		err_n = err_0 - a*q;
		if (m_print_level != 0)
		{
			feLog("%lg (%lg)\n", sqrt(fabs(err_n)), m_tol*sqrt(err_0));
		}

		// check for convergence
		if (sqrt(fabs(err_n)) < m_tol*sqrt(err_0)) break;

		// Update R vectors
		if (n < m_maxiter - 1)
		{
			// update R vectors
			RM = Rm0;
			RC = Rc0;
			for (int i = 0; i <= n; ++i)
			{
				vsubs(RM, Rmu[i], au[i]);
				vsubs(RM, Rmp[i], ap[i]);

				vsubs(RC, Rcu[i], au[i]);
				vsubs(RC, Rcp[i], ap[i]);
			}
		}
	}

	if (m_print_level != 0) feLog("---\n");

	// calculate final solution vector
	yu.assign(Nu, 0.0);
	yp.assign(Np, 0.0);
	for (int i = 0; i<niter; ++i)
	{
		vadds(yu, Yu[i], au[i]);
		vadds(yp, Yp[i], ap[i]);
	}

	// de-normalize the solution
	vscale(yu, Wm);
	vscale(yp, Wc);

	// put it all together
	for (int i = 0; i<Nu; ++i) x[i     ] = yu[i];
	for (int i = 0; i<Np; ++i) x[i + Nu] = yp[i];

	// and ... scene!
	return true;
}

int BIPNSolver::cgsolve(SparseMatrix* K, LinearSolver* PC, std::vector<double>& x, std::vector<double>& b)
{
	m_cg_iters = 0;
	RCICGSolver cg(nullptr);
	if (cg.SetSparseMatrix(K) == false) return -1;

	if (PC) cg.SetLeftPreconditioner(PC);

	cg.SetMaxIterations(m_cg_maxiter);
	cg.SetTolerance(m_cg_tol);
//	cg.SetPrintLevel(1);

	if (cg.PreProcess() == false) return -1;
	if (cg.Factor() == false) return -1;
	cg.BackSolve(&x[0], &b[0]);

	return cg.GetStats().iterations;
}

int BIPNSolver::gmressolve(SparseMatrix* K, LinearSolver* PC, vector<double>& x, vector<double>& b)
{
	FGMRESSolver gmres(nullptr);
	if (gmres.SetSparseMatrix(K) == false) return -1;

	if (PC) gmres.SetLeftPreconditioner(PC);

	gmres.SetMaxIterations(m_gmres_maxiter);
	gmres.SetRelativeResidualTolerance(m_gmres_tol);

	if (gmres.PreProcess() == false) return -1;
	if (gmres.Factor() == false) return -1;
	gmres.BackSolve(&x[0], &b[0]);

	return gmres.GetStats().iterations;
}

#else	// ifdef MKL_ISS
BEGIN_FECORE_CLASS(BIPNSolver, LinearSolver)
	ADD_PARAMETER(m_maxiter, "maxiter");
	ADD_PARAMETER(m_tol, "tol");
	ADD_PARAMETER(m_print_level, "print_level");
	ADD_PARAMETER(m_use_cg, "use_cg");
	ADD_PARAMETER(m_cg_maxiter, "cg_maxiter");
	ADD_PARAMETER(m_cg_tol, "cg_tol");
	ADD_PARAMETER(m_cg_doResidualTest, "cg_check_residual");
	ADD_PARAMETER(m_gmres_maxiter, "gmres_maxiter");
	ADD_PARAMETER(m_gmres_tol, "gmres_tol");
	ADD_PARAMETER(m_gmres_doResidualTest, "gmres_check_residual");
	ADD_PARAMETER(m_gmres_pc, "gmres_precondition");
	ADD_PARAMETER(m_do_jacobi, "do_jacobi");
	ADD_PARAMETER(m_precondition_schur, "precondition_schur");
END_FECORE_CLASS();

BIPNSolver::BIPNSolver(FEModel* fem) : LinearSolver(fem), m_A(0) {}
bool BIPNSolver::PreProcess() { return false; }
bool BIPNSolver::Factor() { return false; }
bool BIPNSolver::BackSolve(double* x, double* b) { return false; }
SparseMatrix* BIPNSolver::CreateSparseMatrix(Matrix_Type ntype) { return 0; }
void BIPNSolver::SetPrintLevel(int n) {}
void BIPNSolver::SetMaxIterations(int n) {}
void BIPNSolver::SetTolerance(double eps) {}
void BIPNSolver::UseConjugateGradient(bool b) {}
void BIPNSolver::SetCGParameters(int maxiter, double tolerance, bool doResidualStoppingTest) {}
void BIPNSolver::SetGMRESParameters(int maxiter, double tolerance, bool doResidualStoppingTest, int precondition) {}
void BIPNSolver::DoJacobiPreconditioner(bool b) {}
void BIPNSolver::SetSchurPreconditioner(int n) {}
bool BIPNSolver::SetSparseMatrix(SparseMatrix* A) { return false; }
#endif
