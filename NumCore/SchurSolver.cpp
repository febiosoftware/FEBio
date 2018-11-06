#include "SchurSolver.h"
#include "FGMRES_ILU0_Solver.h"
#include <FECore/SchurComplement.h>
#include "ILU0_Preconditioner.h"
#include "IncompleteCholesky.h"
#include "HypreGMRESsolver.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEGlobalMatrix.h>

//-----------------------------------------------------------------------------
//! constructor
SchurSolver::SchurSolver(FEModel* fem) : LinearSolver(fem)
{
	m_pA = 0;
	m_tol = 1e-12;
	m_maxiter = 0;
	m_iter = 0;
	m_printLevel = 0;

	m_nsolver = 0;

	m_schurBlock = 0;

	m_PS = nullptr;
	m_schurSolver = nullptr;

	m_buildMassMatrix = false;

	m_bfailMaxIters = true;

	m_bzeroDBlock = false;
}

//-----------------------------------------------------------------------------
//! constructor
SchurSolver::~SchurSolver()
{
}

//-----------------------------------------------------------------------------
void SchurSolver::SetRelativeTolerance(double tol)
{
	m_tol = tol;
}

//-----------------------------------------------------------------------------
// get the iteration count
int SchurSolver::GetIterations() const
{
	return m_iter;
}

//-----------------------------------------------------------------------------
// set the print level
void SchurSolver::SetPrintLevel(int n)
{
	m_printLevel = n;
}

//-----------------------------------------------------------------------------
// set max nr of iterations
void SchurSolver::SetMaxIterations(int n)
{
	m_maxiter = n;
}

//-----------------------------------------------------------------------------
// set convergence tolerance
void SchurSolver::SetConvergenceTolerance(double tol)
{
	m_tol = tol;
}

//-----------------------------------------------------------------------------
void SchurSolver::SetSchurBlock(int n)
{
	m_schurBlock = n;
}

//-----------------------------------------------------------------------------
//! Set the partition
void SchurSolver::SetPartitions(const vector<int>& part)
{
	m_npart = part;
}

//-----------------------------------------------------------------------------
void SchurSolver::UseMassMatrix(bool b)
{
	m_buildMassMatrix = b;
}

//-----------------------------------------------------------------------------
void SchurSolver::SetLinearSolver(int n)
{
	m_nsolver = n;
}

//-----------------------------------------------------------------------------
void SchurSolver::FailOnMaxIterations(bool b)
{
	m_bfailMaxIters = b;
}

//-----------------------------------------------------------------------------
void SchurSolver::ZeroDBlock(bool b)
{
	m_bzeroDBlock = b;
}

//-----------------------------------------------------------------------------
//! Create a sparse matrix
SparseMatrix* SchurSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	if (m_npart.size() != 2) return 0;
	m_pA = new BlockMatrix();
	m_pA->Partition(m_npart, ntype, (m_nsolver == 0 ? 1 : 0));
	return m_pA;
}

//-----------------------------------------------------------------------------
//! set the sparse matrix
bool SchurSolver::SetSparseMatrix(SparseMatrix* A)
{
	m_pA = dynamic_cast<BlockMatrix*>(A);
	if (m_pA == 0) return false;
	return true;
}

//-----------------------------------------------------------------------------
//! Preprocess 
bool SchurSolver::PreProcess()
{
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// get the number of partitions
	// and make sure we have two
	int NP = m_pA->Partitions();
	if (NP != 2) return false;

	// allocate solvers for diagonal blocks
	if (m_nsolver == 0)
	{
		FGMRES_ILU0_Solver* fgmres = new FGMRES_ILU0_Solver(GetFEModel());
		fgmres->SetMaxIterations(m_maxiter);
		fgmres->SetPrintLevel(m_printLevel == 3 ? 0 : m_printLevel);
		fgmres->SetResidualTolerance(m_tol);
		fgmres->FailOnMaxIterations(false);
		m_solver = fgmres;
	}
	else
	{
		HypreGMRESsolver* fgmres = new HypreGMRESsolver(GetFEModel());
		fgmres->SetMaxIterations(m_maxiter);
		fgmres->SetPrintLevel(m_printLevel == 3 ? 0 : m_printLevel);
		fgmres->SetConvergencTolerance(m_tol);
		m_solver = fgmres;
	}

	if (m_schurBlock == 0)
	{
		BlockMatrix::BLOCK& A = m_pA->Block(0, 0);
		m_solver->SetSparseMatrix(A.pA);
	}
	else
	{
		BlockMatrix::BLOCK& D = m_pA->Block(1, 1);
		m_solver->SetSparseMatrix(D.pA);
	}
	if (m_solver->PreProcess() == false) return false;

	FGMRESSolver* fgmres2 = new FGMRESSolver(GetFEModel());
	fgmres2->SetPrintLevel(m_printLevel == 3 ? 2 : m_printLevel);
	if (m_maxiter > 0) fgmres2->SetMaxIterations(m_maxiter);
	fgmres2->SetResidualTolerance(m_tol);
	fgmres2->FailOnMaxIterations(m_bfailMaxIters);

	m_schurSolver = fgmres2;

	m_iter = 0;

	return true;
}

//-----------------------------------------------------------------------------
//! Factor matrix
bool SchurSolver::Factor()
{
	// factor the diagonal matrix
	if (m_solver->Factor() == false) return false;

	if ((m_schurBlock == 0) && (m_buildMassMatrix && (m_PS == nullptr)))
	{
		CompactSymmMatrix* M = new CompactSymmMatrix(1);
		if (BuildMassMatrix(M) == false) return false;

		// We do a LU factorization
		m_PS = new IncompleteCholesky(GetFEModel());
		if (m_PS->Create(M) == false) return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Backsolve the linear system
bool SchurSolver::BackSolve(double* x, double* b)
{
	// get the partition sizes
	int n0 = m_pA->PartitionEquations(0);
	int n1 = m_pA->PartitionEquations(1);

	// Get the blocks
	BlockMatrix::BLOCK& A = m_pA->Block(0, 0);
	BlockMatrix::BLOCK& B = m_pA->Block(0, 1);
	BlockMatrix::BLOCK& C = m_pA->Block(1, 0);
	BlockMatrix::BLOCK& D = m_pA->Block(1, 1);

	// split right hand side in two
	vector<double> F(n0), G(n1);
	for (int i = 0; i<n0; ++i) F[i] = b[i];
	for (int i = 0; i<n1; ++i) G[i] = b[i + n0];

	double normF = l2_norm(F);
	double normG = l2_norm(G);

	// solution vectors
	vector<double> u(n0, 0.0);
	vector<double> v(n1, 0.0);

	bool bconv = false;
	if (m_schurBlock == 0)
	{
		// step 1: solve Ay = F
		vector<double> y(n0);
		if (m_printLevel != 0) fprintf(stderr, "----------------------\nstep 1:\n");
		if (m_solver->BackSolve(y, F) == false) return false;

		// step 2: calculate H = Cy - G
		vector<double> H(n1);
		C.vmult(y, H);
		H -= G;

		// step 3: Solve Sv = H
		SchurComplement S(m_solver, B.pA, C.pA, (m_bzeroDBlock ? nullptr : D.pA));
		if (m_printLevel != 0) fprintf(stderr, "step 3:\n");
		bconv = m_schurSolver->Solve(S, v, H, m_PS);
		if (bconv == false) return false;

		// step 4: calculate L = F - Bv
		vector<double> tmp(n0);
		B.vmult(v, tmp);
		vector<double> L = F - tmp;

		// step 5: solve Au = L
		if (m_printLevel != 0) fprintf(stderr, "step 5:\n");
		m_solver->BackSolve(u, L);
	}
	else
	{
		// step 1: solve Dy = G
		vector<double> y(n1);
		if (m_printLevel != 0) fprintf(stderr, "----------------------\nstep 1:\n");
		if (m_solver->BackSolve(y, G) == false) return false;

		// step 2: calculate H = F - By
		vector<double> H(n0);
		B.vmult(y, H);
		H -= F;

		// step 3: Solve Sv = H
		SchurComplement2 S(A.pA, B.pA, C.pA, m_solver);
		FGMRESSolver fgmres(GetFEModel());
		fgmres.SetPrintLevel(m_printLevel);
		if (m_maxiter > 0) fgmres.SetMaxIterations(m_maxiter);
		fgmres.SetResidualTolerance(m_tol);
		if (m_printLevel != 0) fprintf(stderr, "step 3:\n");
		bconv = fgmres.Solve(&S, u, H);
		if (bconv == false) return false;

		// step 4: calculate L = G - Cu
		vector<double> tmp(n1);
		C.vmult(u, tmp);
		vector<double> L = G - tmp;

		// step 5: solve Dv = L
		if (m_printLevel != 0) fprintf(stderr, "step 5:\n");
		m_solver->BackSolve(v, L);
	}

	// put it back together
	for (int i = 0; i<n0; ++i) x[i] = u[i];
	for (int i = 0; i<n1; ++i) x[i + n0] = v[i];

	return bconv;
}

//-----------------------------------------------------------------------------
//! Clean up
void SchurSolver::Destroy()
{
	m_solver->Destroy();
}


//-----------------------------------------------------------------------------
bool SchurSolver::BuildMassMatrix(CompactSymmMatrix* M)
{
	FEModel* fem = GetFEModel();
	FEMesh& mesh = fem->GetMesh();

	// get number of equations
	int N0 = m_pA->Block(0, 0).Rows();
	int N = m_pA->Block(1, 1).Rows();

	// this is the degree of freedom in the LM arrays that we need
	// TODO: This is hard coded for fluid problems
	const int edofs = 4; // nr of degrees per element node
	const int dof = 3;	// degree of freedom we need

	// build the global matrix
	vector<vector<int> > LM;
	vector<int> lm, lme;
	SparseMatrixProfile MP(N, N);
	MP.CreateDiagonal();
	int ND = mesh.Domains();
	for (int i = 0; i < ND; ++i)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(i));
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j)
		{
			FESolidElement& el = dom.Element(j);
			int neln = el.Nodes();
			// get the equation numbers
			dom.UnpackLM(el, lme);

			// we don't want equation numbers below N0
			lm.resize(neln);
			for (int i = 0; i < neln; ++i)
			{
				lm[i] = lme[i*edofs + dof];
				if (lm[i] >= 0) lm[i] -= N0;
				else if (lm[i] < -1) lm[i] = -lm[i]-2 - N0;
			}

			LM.push_back(lm);
			MP.UpdateProfile(LM, 1);
			LM.clear();
		}
	}
	M->Create(MP);
	M->Zero();

	// build the mass matrix
	matrix me;
	double density = 1.0;
	int n = 0;
	for (int i = 0; i < ND; ++i)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(i));
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j, ++n)
		{
			FESolidElement& el = dom.Element(j);

			// Get the current element's data
			const int nint = el.GaussPoints();
			const int neln = el.Nodes();
			const int ndof = neln;

			// weights at gauss points
			const double *gw = el.GaussWeights();

			matrix me(ndof, ndof);
			me.zero();

			// calculate element stiffness matrix
			for (int n = 0; n<nint; ++n)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(n);
				double Dn = density;

				// shape functions
				double* H = el.H(n);

				// Jacobian
				double J0 = dom.detJ0(el, n)*gw[n];

				for (int i = 0; i<neln; ++i)
					for (int j = 0; j<neln; ++j)
					{
						double mab = Dn*H[i] * H[j] * J0;
						me[i][j] += mab;
					}
			}

			// get the equation numbers
			dom.UnpackLM(el, lme);

			// we don't want equation numbers below N0
			lm.resize(neln);
			for (int i = 0; i < neln; ++i)
			{
				lm[i] = lme[i*edofs + dof];
				if (lm[i] >= 0) lm[i] -= N0;
				else if (lm[i] < -1) lm[i] = -lm[i]-2 - N0;
			}

			M->Assemble(me, lm);
		}
	}

	return true;
}
