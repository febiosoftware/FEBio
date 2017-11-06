#include "stdafx.h"
#include "BIPNSolver.h"
#include <FECore/vector.h>

#ifdef MKL_ISS

// We must undef PARDISO since it is defined as a function in mkl_solver.h
#ifdef PARDISO
#undef PARDISO
#endif

#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"

// constructor
BIPNSolver::BIPNSolver() : m_A(0)
{
	m_print_level = 0;
	m_maxiter = 10;
	m_split = -1;	// invalid! Must be set by user
	m_tol = 1e-6;

	m_cg_maxiter = 0;
	m_cg_tol = 1e-6;

	m_gmres_maxiter = 0;
	m_gmres_tol = 1e-6;

	m_use_cg = true;
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

// set the row/column index that defines the individual blocks
void BIPNSolver::SetPartition(int n)
{
	m_split = n;
}

// Use CG for step 2 or not
void BIPNSolver::UseConjugateGradient(bool b)
{
	m_use_cg = b;
}

//! Return a sparse matrix compatible with this solver
SparseMatrix* BIPNSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	// make sure we can support this matrix
	if (ntype != Matrix_Type::REAL_UNSYMMETRIC) return 0;

	// allocate new matrix
	if (m_A) delete m_A;
	m_A = new CompactUnSymmMatrix(0, true);

	// and return
	return m_A;
}

// allocate storage
bool BIPNSolver::PreProcess()
{
	// make sure we have a matrix
	if (m_A == 0) return false;
	CompactUnSymmMatrix& A = *m_A;

	// get the number of equations
	int N = A.Size();

	// make sure the split is valid
	if ((m_split <= 0) || (m_split >= N-1)) return false;
	int Nu = m_split;
	int Np = N - Nu;

	// allocate pre-conditioner
	m_W.resize(N);
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

	return true;
}

//! Pre-condition the matrix
bool BIPNSolver::Factor()
{
	// make sure we have a matrix
	if (m_A == 0) return false;
	CompactUnSymmMatrix& A = *m_A;

	// get the number of equations
	int N = A.Size();
	int Nu = m_split;
	int Np = N - Nu;

	// calculate preconditioner
	assert(N == (int) m_W.size());
	for (int i=0; i<N; ++i) 
	{
		double di = fabs(A.diag(i));
		assert(di != 0.0);
		m_W[i] = (di != 0.0 ? 1.0 / sqrt(di) : 1.0);
	}

	// split in two
	for (int i=0; i<Nu; ++i) Wm[i] = m_W[i];
	for (int i=0; i<Np; ++i) Wc[i] = m_W[i + Nu];

	// normalize the matrix
	A.scale(m_W, m_W);

	// extract the blocks
	// 
	//     | K | G |
	// A = |---+---|
	//     | D | L |
	//
	A.get( 0,  0, Nu, Nu, K);
	A.get( 0, Nu, Nu, Np, G);
	A.get(Nu,  0, Np, Nu, D);
	A.get(Nu, Nu, Np, Np, L);

	return true;
}

//! Calculate the solution of RHS b and store solution in x
bool BIPNSolver::BackSolve(vector<double>& x, vector<double>& b)
{
	// make sure we have a matrix
	if (m_A == 0) return false;
	CompactUnSymmMatrix& A = *m_A;

	// number of equations
	int N = A.Size();
	int Nu = m_split;
	int Np = N - Nu;

	// normalize RHS
	for (int i = 0; i<Nu; ++i) Rm[i] = Wm[i]*b[i     ];
	for (int i = 0; i<Np; ++i) Rc[i] = Wc[i]*b[i + Nu];

	// initialize
	RM[0] = Rm;
	RC[0] = Rc;

	// calculate initial error
	double err_0 = Rm*Rm + Rc*Rc, err_n = 0.0;

	// do the BIPN iterations 
	int niter = 0;
	for (int n=0; n<m_maxiter; ++n)
	{
		niter++;

		// solve for yu_n (use GMRES): K*yu_n = RM[n]
		gmressolve(yu_n, RM[n], m_gmres_maxiter, m_gmres_tol);

		// compute the corrected residual Rc_n = RC[n] - D*yu_n
		D.multv(yu_n, Rc_n);
		vsub(Rc_n, RC[n], Rc_n);

		// solve for yp_n (use CG): (L + D*G) * yp_n = Rc_n
		if (m_use_cg)
			step2_cgsolve(yp_n, Rc_n, m_cg_maxiter, m_cg_tol);
		else
			step2_gmressolve(yp_n, Rc_n, m_gmres_maxiter, m_gmres_tol);

		// compute corrected residual: Rm_n = RM[n] - G*yp_n
		G.multv(yp_n, Rm_n);
		vsub(Rm_n, RM[n], Rm_n);

		// solve for yu_n (use GMRES): K*yu_n = Rm_n
		gmressolve(yu_n, Rm_n, m_gmres_maxiter, m_gmres_tol);

		// calculate temp vectors
		K.multv(yu_n, Rmu[n]);	// Rmu[n] = K*yu_n;
		G.multv(yp_n, Rmp[n]);	// Rmp[n] = G*yp_n;
		D.multv(yu_n, Rcu[n]);	// Rcu[n] = D*yu_n;
		L.multv(yp_n, Rcp[n]);	// Rcp[n] = L*yp_n;

		// store solution candidates
		Yu[n] = yu_n;
		Yp[n] = yp_n;

		// number of coefficients to determine
		const int m = 2 * (n + 1);

		// setup Q matrix
		matrix Q(m, m);
		vector<double> q(m);

		// fill in the blocks of Q and q
		// TODO: There is an optimization here
		for (int i = 0; i <= n; ++i)
		{
			for (int j = 0; j <= n; ++j)
			{
				Q(i        , j        ) = Rmu[i] * Rmu[j] + Rcu[i] * Rcu[j];
				Q(i        , j + n + 1) = Rmu[i] * Rmp[j] + Rcu[i] * Rcp[j];
				Q(i + n + 1, j        ) = Rmp[i] * Rmu[j] + Rcp[i] * Rcu[j];
				Q(i + n + 1, j + n + 1) = Rmp[i] * Rmp[j] + Rcp[i] * Rcp[j];
			}

			q[i] = Rm*Rmu[i] + Rc*Rcu[i];
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
			printf("BIPN error %d = %lg\n", n, sqrt(fabs(err_n)));
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

bool BIPNSolver::step2_cgsolve(std::vector<double>& x, std::vector<double>& b, int maxiter, double tol)
{
	// make sure we have a matrix
	if (m_A == 0) return false;
	CompactUnSymmMatrix& A = *m_A;

	// make sure this is a symmetric matrix
	int NEQ = A.Size();
	int Nu = m_split;
	int Np = NEQ - Nu;
	int n = Np;

	// output parameters
	MKL_INT rci_request;
	MKL_INT ipar[128];
	double dpar[128];
	double* tmp = &cg_tmp[0];

	// initialize parameters
	dcg_init(&n, &x[0], &b[0], &rci_request, ipar, dpar, tmp);
	if (rci_request != 0) return false;

	// set the desired parameters:
	if (maxiter > 0) ipar[4] = maxiter;		// set max iterations
	ipar[ 8] = 1;		// do residual stopping test
	ipar[ 9] = 0;		// do not do user-defined stopping test
	ipar[10] = 0;		// do or do not precondition
	dpar[ 0] = 1e-5;		// set relative residual tolerance (default 1e-6)

	// check the consistency of the newly set parameters
	dcg_check(&n, &x[0], &b[0], &rci_request, ipar, dpar, tmp);
	if (rci_request != 0) return false;

	// loop until converged
	bool bsuccess = false;
	bool bdone = false;
	do
	{
		// compute the solution by RCI
		dcg(&n, &x[0], &b[0], &rci_request, ipar, dpar, tmp);

//		printf("CG residual : %lg\n", dpar[4]);

		switch (rci_request)
		{
		case 0: // solution converged! 
			bsuccess = true;
			bdone = true;
			break;
		case 1: // compute vector A*tmp[0] and store in tmp[n]
		{
			// In this case, the matrix is defined as S = L + D*G
			// where it is assumed that D = G^t

			// multiply tmp with G and store in du
			G.multv(tmp, &du[0]);

			// multiply du with D and store in tmp+n
			D.multv(&du[0], tmp + n);

			// multiply tmp with L and store in dp
			L.multv(tmp, &dp[0]);

			// now, add dp to tmp + n
			for (int i=0; i<Np; ++i) tmp[n + i] += dp[i];
		}
		break;
		default:
			bsuccess = false;
			bdone = true;
			break;
		}
	}
	while (!bdone);

	// get convergence information
	int niter;
	dcg_get(&n, &x[0], &b[0], &rci_request, ipar, dpar, tmp, &niter);

	if (m_print_level != 0)
	{
		printf("CG iterations: %d\n", niter);
	}

	assert(bsuccess);
	return bsuccess;
}


bool BIPNSolver::step2_gmressolve(vector<double>& x, vector<double>& b, int maxiter, double tol)
{
	// make sure we have a matrix
	if (m_A == 0) return false;
	CompactUnSymmMatrix& A = *m_A;

	// make sure this is a symmetric matrix
	int NEQ = A.Size();
	int Nu = m_split;
	int Np = NEQ - Nu;
	int N = Np;

	// data allocation
	MKL_INT ipar[128];
	double dpar[128];
	MKL_INT RCI_request;
	int M = (N < 150 ? N : 150); // this is the default value of par[15] (i.e. par[14] in C)
	if (maxiter > 0) M = maxiter;

	// allocate temp storage
	double* tmp = &cg_tmp[0];

	MKL_INT ivar = N;

	// initialize the solver
	dfgmres_init(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// Set the desired parameters:
	ipar[4] = M;	// max number of iterations
	ipar[8] = 1;		// do residual stopping test
	ipar[9] = 0;		// do not request for the user defined stopping test
	ipar[11] = 1;		// do the check of the norm of the next generated vector automatically
	//	dpar[0]=1.0E-3;	// set the relative tolerance to 1.0D-3 instead of default value 1.0D-6

	// Check the correctness and consistency of the newly set parameters
	dfgmres_check(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// solve the problem
	bool bdone = false;
	bool bconverged = false;
	while (!bdone)
	{
		// compute the solution via FGMRES
		dfgmres(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, tmp);

		switch (RCI_request)
		{
		case 0: // solution converged. 
			bdone = true;
			bconverged = true;
			break;
		case 1:
		{
			// In this case, the matrix is defined as S = L + D*G
			// where it is assumed that D = G^t

			// multiply tmp with G and store in du
			G.multv(&tmp[ipar[21] - 1], &du[0]);

			// multiply du with D and store in tmp
			D.multv(&du[0], &tmp[ipar[22] - 1]);

			// multiply tmp with L and store in dp
			L.multv(&tmp[ipar[21] - 1], &dp[0]);

			// now, add dp to tmp
			for (int i = 0; i<N; ++i) tmp[ipar[22] -1 + i] += dp[i];
		}
		break;
		default:	// something went wrong
			bdone = true;
			bconverged = false;
		}
	}

	// get the solution. 
	MKL_INT itercount;
	dfgmres_get(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, tmp, &itercount);

	assert(bconverged);
	return bconverged;
}

bool BIPNSolver::gmressolve(vector<double>& x, vector<double>& b, int maxiter, double tol)
{
	// make sure we have a matrix
	if (m_A == 0) return false;
	CompactUnSymmMatrix& A = *m_A;

	// make sure this is a symmetric matrix
	int NEQ = A.Size();
	int Nu = m_split;
	int Np = NEQ - Nu;
	int N = Nu;

	// data allocation
	MKL_INT ipar[128];
	double dpar[128];
	MKL_INT RCI_request;
	int M = (N < 150 ? N : 150); // this is the default value of par[15] (i.e. par[14] in C)
	if (maxiter > 0) M = maxiter;

	// allocate temp storage
	double* tmp = &gmres_tmp[0];

	MKL_INT ivar = N;

	// initialize the solver
	dfgmres_init(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// Set the desired parameters:
	ipar[4] = M;	// max number of iterations
	ipar[8] = 1;		// do residual stopping test
	ipar[9] = 0;		// do not request for the user defined stopping test
	ipar[11] = 1;		// do the check of the norm of the next generated vector automatically
	//	dpar[0]=1.0E-3;	// set the relative tolerance to 1.0D-3 instead of default value 1.0D-6

	// Check the correctness and consistency of the newly set parameters
	dfgmres_check(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// solve the problem
	bool bdone = false;
	bool bconverged = false;
	while (!bdone)
	{
		// compute the solution via FGMRES
		dfgmres(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, tmp);

		switch (RCI_request)
		{
		case 0: // solution converged. 
			bdone = true;
			bconverged = true;
			break;
		case 1:
		{
			// multiply with matrix
			K.multv(&tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);
		}
		break;
		default:	// something went wrong
			bdone = true;
			bconverged = false;
		}
	}

	// get the solution. 
	MKL_INT itercount;
	dfgmres_get(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, tmp, &itercount);

	assert(bconverged);
	return bconverged;
}

#else	// ifdef MKL_ISS

BIPNSolver::BIPNSolver() : m_K(0) {}
bool BIPNSolver::Factor() { return false; }
bool BIPNSolver::BackSolve(vector<double>& x, vector<double>& b) { return false; }
SparseMatrix* BIPNSolver::CreateSparseMatrix(Matrix_Type ntype) { return 0; }

#endif
