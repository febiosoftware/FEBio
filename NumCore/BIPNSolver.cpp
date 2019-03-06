#include "stdafx.h"
#include "BIPNSolver.h"
#include <FECore/vector.h>
#include <FECore/FEModel.h>
#include "RCICGSolver.h"
#include "FGMRES_ILU0_Solver.h"
#include <FECore/SchurComplement.h>
#include <FECore/log.h>
#include "SchurSolver.h" // for PCSolver
#include "IncompleteCholesky.h"

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
	m_gmres_ilu0 = false;

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
void BIPNSolver::SetGMRESParameters(int maxiter, double tolerance, bool doResidualStoppingTest, bool precondition)
{
	m_gmres_maxiter        = maxiter;
	m_gmres_tol            = tolerance;
	m_gmres_doResidualTest = doResidualStoppingTest;
	m_gmres_ilu0           = precondition;
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

	// allocate new matrix
	if (m_A) delete m_A;
	m_A = new BlockMatrix();
	m_A->Partition(m_part, ntype);

	// and return
	return m_A;
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

	Rm.resize(Nu);
	Rc.resize(Np);

	Rm_n.resize(Nu);
	Rc_n.resize(Np);

	Yu.resize(m_maxiter, std::vector<double>(Nu));
	Yp.resize(m_maxiter, std::vector<double>(Np));

	RM.resize(m_maxiter, std::vector<double>(Nu));
	RC.resize(m_maxiter, std::vector<double>(Np));

	Rmu.resize(m_maxiter, std::vector<double>(Nu));
	Rmp.resize(m_maxiter, std::vector<double>(Nu));
	Rcu.resize(m_maxiter, std::vector<double>(Np));
	Rcp.resize(m_maxiter, std::vector<double>(Np));

	au.resize(Nu);
	ap.resize(Np);

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

	// initialize solver for A block
	if (m_gmres_ilu0)
		m_Asolver = new FGMRES_ILU0_Solver(nullptr);
	else
		m_Asolver = new FGMRESSolver(nullptr);

	m_Asolver->SetMaxIterations(m_gmres_maxiter);
	m_Asolver->SetRelativeResidualTolerance(m_gmres_tol);
	if (m_Asolver->SetSparseMatrix(A.Block(0, 0).pA) == false) return false;

	if (m_Asolver->GetPreconditioner())
	{
		m_Asolver->GetPreconditioner()->SetSparseMatrix(A.Block(0, 0).pA);
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
			PS->SetSparseMatrix(M);
			if (PS->Create() == false) return false;

			m_PS = PS;
		}
		else
		{
			// mass matrix
			CompactSymmMatrix* M = new CompactSymmMatrix(1);
			if (BuildMassMatrix(GetFEModel(), m_A, M, 1.0) == false) return false;

			// We do an incomplete cholesky factorization
			IncompleteCholesky* PS = new IncompleteCholesky(GetFEModel());
			PS->SetSparseMatrix(M);
			if (PS->Create() == false) return false;

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
			double ki = K.pA->diag(i);
			if (ki <= 0.0) return false;
			Wm[i] = 1.0 / sqrt(ki);
		}
		if (L.pA)
		{
			for (int i = 0; i < Np; ++i)
			{
				double li = L.pA->diag(i);
				if (li < 0.0) return false;
				Wc[i] = (li > 0.0 ? 1.0 / sqrt(li) : 1.0);
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
	for (int i = 0; i<Nu; ++i) Rm[i] = Wm[i]*b[i     ];
	for (int i = 0; i<Np; ++i) Rc[i] = Wc[i]*b[i + Nu];

	// initialize
	RM[0] = Rm;
	RC[0] = Rc;

	// calculate initial error
	double err_0 = Rm*Rm + Rc*Rc, err_n = 0.0;

	// setup the Schur complement
	PCSolver PC(nullptr);
	DiagonalPreconditioner DPC(nullptr);
	DPC.SetSparseMatrix(K.pA);
	DPC.Create();
	PC.SetPreconditioner(&DPC);
	SchurComplement S(&PC, G.pA, D.pA, L.pA);

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
		m_Asolver->BackSolve(&yu_n[0], &(RM[n][0]));
		m_gmres1_iters = m_Asolver->GetStats().iterations;
		if (m_print_level != 0) feLog("%d, ", m_gmres1_iters);

		// compute the corrected residual Rc_n = RC[n] - D*yu_n
		D.vmult(yu_n, Rc_n);
		vsub(Rc_n, RC[n], Rc_n);

		// solve for yp_n (use CG): (L + D*G) * yp_n = Rc_n
		if (m_use_cg)
			m_cg_iters = cgsolve(&S, m_PS, yp_n, Rc_n);
		else
			m_cg_iters = gmressolve(&S, m_PS, yp_n, Rc_n);
		if (m_print_level != 0) feLog("%d, ", m_cg_iters);

		// compute corrected residual: Rm_n = RM[n] - G*yp_n
		G.vmult(yp_n, Rm_n);
		vsub(Rm_n, RM[n], Rm_n);

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

			q[i        ] = Rm*Rmu[i] + Rc*Rcu[i];
			q[i + n + 1] = Rm*Rmp[i] + Rc*Rcp[i];
		}

		// solve for the coefficients
		vector<double> a(m);
		Q.solve(q, a);
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
			RM[n + 1] = Rm;
			RC[n + 1] = Rc;
			for (int i = 0; i <= n; ++i)
			{
				vsubs(RM[n + 1], Rmu[i], au[i]);
				vsubs(RM[n + 1], Rmp[i], ap[i]);

				vsubs(RC[n + 1], Rcu[i], au[i]);
				vsubs(RC[n + 1], Rcp[i], ap[i]);
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

int BIPNSolver::cgsolve(SparseMatrix* K, Preconditioner* PC, std::vector<double>& x, std::vector<double>& b)
{
	m_cg_iters = 0;
	RCICGSolver cg(nullptr);
	if (cg.SetSparseMatrix(K) == false) return -1;

	if (PC) cg.SetPreconditioner(PC);

	cg.SetMaxIterations(m_cg_maxiter);
	cg.SetTolerance(m_cg_tol);
//	cg.SetPrintLevel(1);

	if (cg.PreProcess() == false) return -1;
	if (cg.Factor() == false) return -1;
	cg.BackSolve(&x[0], &b[0]);

	return cg.GetStats().iterations;
}

int BIPNSolver::gmressolve(SparseMatrix* K, Preconditioner* PC, vector<double>& x, vector<double>& b)
{
	FGMRESSolver gmres(nullptr);
	if (gmres.SetSparseMatrix(K) == false) return -1;

	if (PC) gmres.SetPreconditioner(PC);

	gmres.SetMaxIterations(m_gmres_maxiter);
	gmres.SetRelativeResidualTolerance(m_gmres_tol);

	if (gmres.PreProcess() == false) return -1;
	if (gmres.Factor() == false) return -1;
	gmres.BackSolve(&x[0], &b[0]);

	return gmres.GetStats().iterations;
}

#else	// ifdef MKL_ISS

BIPNSolver::BIPNSolver() : LinearSolver(fem), m_A(0) {}
bool BIPNSolver::PreProcess() { return false; }
bool BIPNSolver::Factor() { return false; }
bool BIPNSolver::BackSolve(double* x, double* b) { return false; }
SparseMatrix* BIPNSolver::CreateSparseMatrix(Matrix_Type ntype) { return 0; }
void BIPNSolver::SetPrintLevel(int n) {}
void BIPNSolver::SetMaxIterations(int n) {}
void BIPNSolver::SetPartition(int n) {}
void BIPNSolver::SetTolerance(double eps) {}
void BIPNSolver::UseConjugateGradient(bool b) {}
void BIPNSolver::SetCGParameters(int maxiter, double tolerance, bool doResidualStoppingTest) {}
void BIPNSolver::SetGMRESParameters(int maxiter, double tolerance, bool doResidualStoppingTest, bool precondition) {}

#endif
