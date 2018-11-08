#include "SchurSolver.h"
#include "FGMRES_ILU0_Solver.h"
#include <FECore/SchurComplement.h>
#include "ILU0_Preconditioner.h"
#include "IncompleteCholesky.h"
#include "HypreGMRESsolver.h"
#include "RCICGSolver.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEGlobalMatrix.h>
#include "PardisoSolver.h"

//-----------------------------------------------------------------------------
//! constructor
SchurSolver::SchurSolver(FEModel* fem) : LinearSolver(fem)
{
	m_pA = 0;
	m_reltol = 1e-8;
	m_abstol = 0.0;
	m_maxiter = 0;
	m_iter = 0;
	m_printLevel = 0;

	m_nsolver = 0;
	m_nschurSolver = 0;

	m_PS = nullptr;
	m_schurSolver = nullptr;

	m_bfailMaxIters = true;

	m_bzeroDBlock = false;
}

//-----------------------------------------------------------------------------
//! constructor
SchurSolver::~SchurSolver()
{
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
void SchurSolver::SetRelativeResidualTolerance(double tol)
{
	m_reltol = tol;
}

//-----------------------------------------------------------------------------
void SchurSolver::SetAbsoluteResidualTolerance(double tol)
{
	m_abstol = tol;
}

//-----------------------------------------------------------------------------
//! Set the partition
void SchurSolver::SetPartitions(const vector<int>& part)
{
	m_npart = part;
}

//-----------------------------------------------------------------------------
void SchurSolver::SetLinearSolver(int n)
{
	m_nsolver = n;
}

//-----------------------------------------------------------------------------
void SchurSolver::SetSchurSolver(int n)
{
	m_nschurSolver = n;
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
	m_pA->Partition(m_npart, ntype, (m_nsolver == 2 ? 0 : 1));
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

	// Get the blocks
	BlockMatrix::BLOCK& A = m_pA->Block(0, 0);
	BlockMatrix::BLOCK& B = m_pA->Block(0, 1);
	BlockMatrix::BLOCK& C = m_pA->Block(1, 0);
	BlockMatrix::BLOCK& D = m_pA->Block(1, 1);

	// get the number of partitions
	// and make sure we have two
	int NP = m_pA->Partitions();
	if (NP != 2) return false;

	// build solver for upper diagonal block
	if (m_nsolver == 0)
	{
		FGMRES_ILU0_Solver* fgmres = new FGMRES_ILU0_Solver(GetFEModel());
		fgmres->SetMaxIterations(m_maxiter);
		fgmres->SetPrintLevel(m_printLevel == 3 ? 0 : m_printLevel);
		fgmres->SetRelativeResidualTolerance(m_reltol);
		fgmres->FailOnMaxIterations(false);
		m_solver = fgmres;
	}
	else if (m_nsolver == 1)
	{
		m_solver = new PardisoSolver(GetFEModel());
//		m_solver = new ILU0_Solver(GetFEModel());
	}
	else
	{
		HypreGMRESsolver* fgmres = new HypreGMRESsolver(GetFEModel());
		fgmres->SetMaxIterations(m_maxiter);
		fgmres->SetPrintLevel(m_printLevel == 3 ? 0 : m_printLevel);
		fgmres->SetConvergencTolerance(m_reltol);
		m_solver = fgmres;
	}

	m_solver->SetSparseMatrix(A.pA);
	if (m_solver->PreProcess() == false) return false;

	if ((m_nschurSolver == 0) || (m_nschurSolver == 1))
	{
		// build solver for Schur complement
		FGMRESSolver* fgmres2 = new FGMRESSolver(GetFEModel());
		fgmres2->SetPrintLevel(m_printLevel == 3 ? 2 : m_printLevel);
		if (m_maxiter > 0) fgmres2->SetMaxIterations(m_maxiter);
		fgmres2->SetRelativeResidualTolerance(m_reltol);
		fgmres2->SetAbsoluteResidualTolerance(m_abstol);
		fgmres2->FailOnMaxIterations(m_bfailMaxIters);

		SchurComplement* S = new SchurComplement(m_solver, B.pA, C.pA, (m_bzeroDBlock ? nullptr : D.pA));
		if (fgmres2->SetSparseMatrix(S) == false) { delete S;  return false; }

		// for schurSolver option 1 we add the mass matrix as preconditioner
		if (m_nschurSolver == 1)
		{
			if (m_PS == nullptr)
			{
				CompactSymmMatrix* M = new CompactSymmMatrix(1);
				if (BuildMassMatrix(M) == false) return false;

				// We do a LU factorization
				m_PS = new IncompleteCholesky(GetFEModel());
				if (m_PS->Create(M) == false) return false;
			}
			fgmres2->SetPreconditioner(m_PS);
		}
		m_schurSolver = fgmres2;
	}
	else
	{
		RCICG_ICHOL_Solver* cg = new RCICG_ICHOL_Solver(GetFEModel());
		cg->SetPrintLevel(m_printLevel == 3 ? 2 : m_printLevel);
		if (m_maxiter > 0) cg->SetMaxIterations(m_maxiter);
		cg->SetTolerance(m_reltol);

		CompactSymmMatrix* M = new CompactSymmMatrix(1);
		if (BuildMassMatrix(M, 10.0) == false) return false;

		if (cg->SetSparseMatrix(M) == false) return false;
		m_schurSolver = cg;
	}
	if (m_schurSolver->PreProcess() == false) return false;

	// reset iteration counter
	m_iter = 0;

	return true;
}

//-----------------------------------------------------------------------------
//! Factor matrix
bool SchurSolver::Factor()
{
	if (m_solver->Factor() == false) return false;
	if (m_schurSolver->Factor() == false) return false;
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
	// step 1: solve Ay = F
	vector<double> y(n0);
	if (m_printLevel != 0) fprintf(stderr, "----------------------\nstep 1:\n");
	if (m_solver->BackSolve(y, F) == false) return false;

	// step 2: calculate H = Cy - G
	vector<double> H(n1);
	C.vmult(y, H);
	H -= G;

	// step 3: Solve Sv = H
	if (m_printLevel != 0) fprintf(stderr, "step 3:\n");
	bconv = m_schurSolver->BackSolve(v, H);
	if (bconv == false) return false;

	// step 4: calculate L = F - Bv
	vector<double> tmp(n0);
	B.vmult(v, tmp);
	vector<double> L = F - tmp;

	// step 5: solve Au = L
	if (m_printLevel != 0) fprintf(stderr, "step 5:\n");
	m_solver->BackSolve(u, L);

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
bool SchurSolver::BuildMassMatrix(CompactSymmMatrix* M, double scale)
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
						me[i][j] += mab*scale;
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
