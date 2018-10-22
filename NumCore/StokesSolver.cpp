#include "StokesSolver.h"
#include "RCICGSolver.h"
#include "SchurComplement.h"
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>

//-----------------------------------------------------------------------------
//! constructor
StokesSolver::StokesSolver(FEModel* fem) : LinearSolver(fem)
{
	m_pA = nullptr;
	m_tol = 1e-12;
	m_maxiter = 0;
	m_iter = 0;
	m_printLevel = 0;
	m_PA = nullptr;
	m_PS = nullptr;
	m_buildMassMatrix = false;
}

//-----------------------------------------------------------------------------
//! constructor
StokesSolver::~StokesSolver()
{
}

//-----------------------------------------------------------------------------
void StokesSolver::SetRelativeTolerance(double tol)
{
	m_tol = tol;
}

//-----------------------------------------------------------------------------
// get the iteration count
int StokesSolver::GetIterations() const
{
	return m_iter;
}

//-----------------------------------------------------------------------------
// set the print level
void StokesSolver::SetPrintLevel(int n)
{
	m_printLevel = n;
}

//-----------------------------------------------------------------------------
// set max nr of iterations
void StokesSolver::SetMaxIterations(int n)
{
	m_maxiter = n;
}

//-----------------------------------------------------------------------------
// set convergence tolerance
void StokesSolver::SetConvergenceTolerance(double tol)
{
	m_tol = tol;
}

//-----------------------------------------------------------------------------
//! Set the partition
void StokesSolver::SetPartitions(const vector<int>& part)
{
	m_npart = part;
}

//-----------------------------------------------------------------------------
//! Create a sparse matrix
SparseMatrix* StokesSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	if (ntype != REAL_SYMMETRIC) return 0;

	if (m_npart.size() != 2) return 0;
	m_pA = new BlockMatrix();
	m_pA->Partition(m_npart, ntype);
	return m_pA;
}

//-----------------------------------------------------------------------------
//! set the sparse matrix
bool StokesSolver::SetSparseMatrix(SparseMatrix* A)
{
	m_pA = dynamic_cast<BlockMatrix*>(A);
	if (m_pA == 0) return false;
	return true;
}

//-----------------------------------------------------------------------------
//! Preprocess 
bool StokesSolver::PreProcess()
{
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// get the number of partitions
	// and make sure we have two
	int NP = m_pA->Partitions();
	if (NP != 2) return false;

	// allocate solvers for diagonal blocks
	RCICGSolver* cg = new RCICGSolver(GetFEModel());
	cg->SetPreconditioner(m_PA = new DiagonalPreconditioner);
	cg->SetMaxIterations(m_maxiter);
	cg->SetPrintLevel(m_printLevel);
	m_solver = cg;
	BlockMatrix::BLOCK& Bi = m_pA->Block(0, 0);
	m_solver->SetSparseMatrix(Bi.pA);
	if (m_solver->PreProcess() == false) return false;

	m_iter = 0;

	return true;
}

//-----------------------------------------------------------------------------
//! Factor matrix
bool StokesSolver::Factor()
{
	if (m_solver->Factor() == false) return false;

	// factor the diagonal matrix
	if (m_PA)
	{
		BlockMatrix::BLOCK& Bi = m_pA->Block(0, 0);
		m_PA->Create(Bi.pA);
	}

	// build the mass matrix
	Preconditioner* Minv = 0;
	if (m_buildMassMatrix && (m_PS == nullptr))
	{
		CompactSymmMatrix* M = new CompactSymmMatrix;
		if (BuildMassMatrix(M) == false) return false;

		// We do a LU factorization
		m_PS = new LUPreconditioner;
		m_PS->Create(M);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Backsolve the linear system
bool StokesSolver::BackSolve(double* x, double* b)
{
	// get the partition sizes
	int n0 = m_pA->PartitionEquations(0);
	int n1 = m_pA->PartitionEquations(1);

	// Get the blocks
	BlockMatrix::BLOCK& A = m_pA->Block(0, 0);
	BlockMatrix::BLOCK& B = m_pA->Block(0, 1);
	BlockMatrix::BLOCK& C = m_pA->Block(1, 0);

	// split right hand side in two
	vector<double> F(n0), G(n1);
	for (int i=0; i<n0; ++i) F[i] = b[i];
	for (int i=0; i<n1; ++i) G[i] = b[i + n0];

	// step 1: solve Ay = F
	vector<double> y(n0);
	if (m_printLevel == 2) fprintf(stdout, "----------------------\nstep 1:\n");
	if (m_solver->BackSolve(y, F) == false) return false;

	// step 2: calculate H = Cy - G
	vector<double> H(n1);
	C.vmult(y, H);
	H -= G;

	// step 3: Solve Sv = H
	SchurComplement S(m_solver, B.pA, C.pA);
	vector<double> v(n1);
	RCICGSolver cg(GetFEModel());
	cg.SetPrintLevel(m_printLevel);
	if (m_maxiter > 0) cg.SetMaxIterations(m_maxiter);
	if (m_printLevel == 2) fprintf(stdout, "step 3:\n");
	if (cg.Solve(&S, v, H, m_PS) == false) return false;

	// step 4: calculate L = F - Bv
	vector<double> tmp(n0);
	B.vmult(v, tmp);
	vector<double> L = F - tmp;

	// step 5: solve Au = L
	vector<double> u(n0, 0.0);
	if (m_printLevel == 2) fprintf(stdout, "step 5:\n");
	if (m_solver->BackSolve(u, L) == false) return false;

	// put it back together
	for (int i=0; i<n0; ++i) x[i   ] = u[i];
	for (int i=0; i<n1; ++i) x[i+n0] = v[i];

	return true;
}

//-----------------------------------------------------------------------------
//! Clean up
void StokesSolver::Destroy()
{
	m_solver->Destroy();
}

//-----------------------------------------------------------------------------
bool StokesSolver::BuildMassMatrix(CompactSymmMatrix* M)
{
	FEModel* fem = GetFEModel();
	FEMesh& mesh = fem->GetMesh();
	FEGlobalMatrix G(M);

	// get number of equations
	int N = m_pA->Rows();

	// build the global matrix
	if (G.Create(mesh, N) == false) return false;

	// zero it initially
	G.Zero();

	// build the mass matrix
	matrix me;
	double density = 1.0;
	int ND = mesh.Domains();
	for (int i = 0; i < ND; ++i)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(i));
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j)
		{
			FESolidElement& el = dom.Element(j);
			int neln = el.Nodes();

			int ndof = 3 * neln;
			me.resize(ndof, ndof);
			me.zero();

			double* w = el.GaussWeights();
			int nint = el.GaussPoints();
			for (int n = 0; n < nint; ++n)
			{
				double* H = el.H(n);
				double detJ = dom.detJ0(el, n);

				for (int a=0; a<neln; ++a)
					for (int b = 0; b < neln; ++b)
					{
						double kab = (H[a] * H[b])*(density * detJ * w[n]);

						me[3*a    ][3*b    ] += kab;
						me[3*a + 1][3*b + 1] += kab;
						me[3*a + 2][3*b + 2] += kab;
					}
			}

			vector<int> lm(ndof);
			for (int n=0; n<neln; ++n)
			{ 
				FENode& node = mesh.Node(el.m_node[n]);
				lm[3*n    ] = node.m_ID[0];
				lm[3*n + 1] = node.m_ID[1];
				lm[3*n + 2] = node.m_ID[2];
			}

			G.Assemble(me, lm);
		}
	}

	return true;
}
