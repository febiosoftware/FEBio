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
#include "SchurSolver.h"
#include "ILU0_Preconditioner.h"
#include <FECore/SchurComplement.h>
#include "ILU0_Preconditioner.h"
#include "IncompleteCholesky.h"
#include "HypreGMRESsolver.h"
#include "RCICGSolver.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEGlobalMatrix.h>
#include "PardisoSolver.h"
#include "BoomerAMGSolver.h"
#include "FGMRESSolver.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
bool BuildDiagonalMassMatrix(FEModel* fem, BlockMatrix* K, CompactSymmMatrix* M, double scale)
{
	FEMesh& mesh = fem->GetMesh();

	// get number of equations
	int N0 = K->Block(0, 0).Rows();
	int N = K->Block(1, 1).Rows();

	// build the global matrix
	SparseMatrixProfile MP(N, N);
	MP.CreateDiagonal();
	M->Create(MP);
	M->Zero();

	// build the mass matrix
	matrix me;
	double density = 1.0;
	int n = 0;
	for (int i = 0; i < mesh.Domains(); ++i)
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
			double Me = 0.0;
			for (int n = 0; n<nint; ++n)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(n);
				double Dn = density;

				// Jacobian
				double Jw = dom.detJ0(el, n)*gw[n];

				Me += Dn*Jw;
			}

			int lm = el.m_lm - N0;
			M->set(lm, lm, Me);
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool BuildMassMatrix(FEModel* fem, BlockMatrix* K, CompactSymmMatrix* M, double scale)
{
	FEMesh& mesh = fem->GetMesh();

	// get number of equations
	int N0 = K->Block(0, 0).Rows();
	int N = K->Block(1, 1).Rows();

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
				else if (lm[i] < -1) lm[i] = -lm[i] - 2 - N0;
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
				else if (lm[i] < -1) lm[i] = -lm[i] - 2 - N0;
			}

			M->Assemble(me, lm);
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(SchurSolver, LinearSolver)
	ADD_PARAMETER(m_printLevel   , "print_level");
	ADD_PARAMETER(m_schurBlock   , "schur_block");
	ADD_PARAMETER(m_doJacobi     , "do_jacobi");
	ADD_PARAMETER(m_bzeroDBlock  , "zero_D_block");
	ADD_PARAMETER(m_nSchurPreC   , "schur_pc");

	ADD_PROPERTY(m_Asolver, "A_solver");
	ADD_PROPERTY(m_schurSolver, "schur_solver");
	ADD_PROPERTY(m_SchurAsolver, "schur_A_solver");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
SchurSolver::SchurSolver(FEModel* fem) : LinearSolver(fem)
{
	// default parameters
	m_printLevel = 0;
	m_bzeroDBlock = false;
	m_schurBlock = 0;		// Use Schur complement S\A

	// default solution strategy
	m_nSchurPreC = 0;		// no preconditioner
	m_iter = 0;
	m_doJacobi = false;

	// initialize pointers
	m_pK = nullptr;
	m_Asolver = nullptr;
	m_PS = nullptr;
	m_schurSolver = nullptr;
	m_SchurAsolver = nullptr;
	m_Acopy = nullptr;
}

//-----------------------------------------------------------------------------
//! constructor
SchurSolver::~SchurSolver()
{
	if (m_Acopy) delete m_Acopy;
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
void SchurSolver::SetSchurPreconditioner(int n)
{
	m_nSchurPreC = n;
}

//-----------------------------------------------------------------------------
void SchurSolver::ZeroDBlock(bool b)
{
	m_bzeroDBlock = b;
}

//-----------------------------------------------------------------------------
// Sets which Schur complement to use
// n = 0: use S\A
// n = 1: use S\D
void SchurSolver::SetSchurBlock(int n)
{
	assert((n == 0) || (n == 1));
	m_schurBlock = n;
}

//-----------------------------------------------------------------------------
void SchurSolver::DoJacobiPreconditioning(bool b)
{
	m_doJacobi = b;
}

//-----------------------------------------------------------------------------
//! Create a sparse matrix
SparseMatrix* SchurSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	if (m_part.size() != 2) return 0;
	m_pK = new BlockMatrix();
	m_pK->Partition(m_part, ntype, 1);

	// we want the A solver to define the blocks
	SparseMatrix* K11 = m_Asolver->CreateSparseMatrix(ntype);
	if (K11) 
	{
		CompactMatrix* A = dynamic_cast<CompactMatrix*>(K11);
		if (A == nullptr) { delete K11; delete m_pK; return nullptr; }
		delete m_pK->Block(0, 0).pA; m_pK->Block(0, 0).pA = A; 
	}
	else
	{
		delete m_pK;
		m_pK = nullptr;
	}

	return m_pK;
}

//-----------------------------------------------------------------------------
//! set the sparse matrix
bool SchurSolver::SetSparseMatrix(SparseMatrix* A)
{
	m_pK = dynamic_cast<BlockMatrix*>(A);
	if (m_pK == 0) return false;
	return true;
}

//-----------------------------------------------------------------------------
// allocate Schur complement solver
LinearSolver* SchurSolver::BuildSchurPreconditioner(int nopt)
{
	switch (nopt)
	{
	case Schur_PC_NONE:
		// no preconditioner selected
		return nullptr;
		break;
	case Schur_PC_DIAGONAL_MASS:
	{
		// diagonal mass matrix
		CompactSymmMatrix* M = new CompactSymmMatrix(1);
		if (BuildDiagonalMassMatrix(GetFEModel(), m_pK, M, 1.0) == false) return nullptr;

		DiagonalPreconditioner* PS = new DiagonalPreconditioner(GetFEModel());
		if (PS->Create(M) == false) return nullptr;

		return PS;
	}
	break;
	case Schur_PC_ICHOL_MASS:
	{
		// mass matrix
		CompactSymmMatrix* M = new CompactSymmMatrix(1);
		if (BuildMassMatrix(GetFEModel(), m_pK, M, 1.0) == false) return nullptr;

		// We do an incomplete cholesky factorization
		IncompleteCholesky* PS = new IncompleteCholesky(GetFEModel());
		if (PS->Create(M) == false) return nullptr;

		return PS;
	}
	break;
	default:
		assert(false);
	};

	return nullptr;
}

//-----------------------------------------------------------------------------
//! Preprocess 
bool SchurSolver::PreProcess()
{
	// make sure we have a matrix
	if (m_pK == 0) return false;

	// Get the blocks
	BlockMatrix::BLOCK& A = m_pK->Block(0, 0);
	BlockMatrix::BLOCK& B = m_pK->Block(0, 1);
	BlockMatrix::BLOCK& C = m_pK->Block(1, 0);
	BlockMatrix::BLOCK& D = m_pK->Block(1, 1);

	// get the number of partitions
	// and make sure we have two
	int NP = m_pK->Partitions();
	if (NP != 2) return false;

	// check the A solver
	if (m_Asolver == nullptr) return false;
	m_Asolver->SetFEModel(GetFEModel());

	// check the schur solver
	if (m_schurSolver == nullptr) return false;
	m_schurSolver->SetFEModel(GetFEModel());

	// build solver for A block in Schur solver
	if (m_SchurAsolver == nullptr) m_SchurAsolver = m_Asolver;
	else m_SchurAsolver->SetFEModel(GetFEModel());

	if (m_schurBlock == 0)
	{
		// Use the Schur complement of A
		m_Asolver->SetSparseMatrix(A.pA);
		if (m_SchurAsolver != m_Asolver)
		{
			SparseMatrix* AforSchur = A.pA;
			m_SchurAsolver->SetSparseMatrix(AforSchur);
		}
		SchurComplementA* S_A = new SchurComplementA(m_SchurAsolver, B.pA, C.pA, (m_bzeroDBlock ? nullptr : D.pA));
		if (m_schurSolver->SetSparseMatrix(S_A) == false) { delete S_A;  return false; }
	}
	else
	{
		// Use the Schur complement of D
		m_Asolver->SetSparseMatrix(D.pA);
		if (m_SchurAsolver != m_Asolver)
		{
			SparseMatrix* DforSchur = D.pA;
			m_SchurAsolver->SetSparseMatrix(DforSchur);
		}
		SchurComplementD* S_D = new SchurComplementD(A.pA, B.pA, C.pA, m_SchurAsolver);
		if (m_schurSolver->SetSparseMatrix(S_D) == false) { delete S_D;  return false; }
	}

	if (m_SchurAsolver != m_Asolver)
	{
		if (m_SchurAsolver->PreProcess() == false) return false;
	}

	if (m_Asolver->PreProcess() == false) return false;

	// build a preconditioner for the schur complement solver
	m_PS = BuildSchurPreconditioner(m_nSchurPreC);
	if (m_PS) m_schurSolver->SetLeftPreconditioner(m_PS);

	if (m_schurSolver->PreProcess() == false) return false;

	// reset iteration counter
	m_iter = 0;

	return true;
}

//-----------------------------------------------------------------------------
//! Factor matrix
bool SchurSolver::Factor()
{
	// Get the blocks
	BlockMatrix::BLOCK& A = m_pK->Block(0, 0);
	BlockMatrix::BLOCK& B = m_pK->Block(0, 1);
	BlockMatrix::BLOCK& C = m_pK->Block(1, 0);
	BlockMatrix::BLOCK& D = m_pK->Block(1, 1);

	// Get the block matrices
	CRSSparseMatrix* MA = dynamic_cast<CRSSparseMatrix*>(A.pA);
	CRSSparseMatrix* MB = dynamic_cast<CRSSparseMatrix*>(B.pA);
	CRSSparseMatrix* MC = dynamic_cast<CRSSparseMatrix*>(C.pA);
	CRSSparseMatrix* MD = dynamic_cast<CRSSparseMatrix*>(D.pA);

	// get the number of partitions
	// and make sure we have two
	int NP = m_pK->Partitions();
	if (NP != 2) return false;

	// get the partition sizes
	int n0 = m_pK->PartitionEquations(0);
	int n1 = m_pK->PartitionEquations(1);

	// calculate the diagonal
	m_Wu.assign(n0, 1.0);
	m_Wp.assign(n1, 1.0);

	if (m_doJacobi)
	{
		for (int i = 0; i < n0; ++i)
		{
			double dii = fabs(MA->diag(i));
//			if (dii < 0) return false;
			if (dii != 0.0) m_Wu[i] = 1.0 / sqrt(dii);
		}
		if (MD)
		{
			for (int i = 0; i < n1; ++i)
			{
				double dii = fabs(MD->diag(i));
//				if (dii < 0) return false;
				if (dii != 0.0) m_Wp[i] = 1.0 / sqrt(dii);
			}
		}

		// scale the matrices
		MA->scale(m_Wu, m_Wu);
		MB->scale(m_Wu, m_Wp);
		MC->scale(m_Wp, m_Wu);
		if (MD) MD->scale(m_Wp, m_Wp);
	}

	// See if we need to copy the A (or D) block
	if (m_Acopy)
	{
		if (m_schurBlock == 0) m_Acopy->CopyValues(A.pA);
		else m_Acopy->CopyValues(D.pA);
	}

	// factor the A block solver
	if (m_Asolver->Factor() == false) return false;

	// factor the A block solver of the Schur complement (if necessary)
	if (m_SchurAsolver && (m_SchurAsolver != m_Asolver))
	{
		if (m_SchurAsolver->Factor() == false) return false;
	}

	// factor the schur complement solver
	if (m_schurSolver->Factor() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! Backsolve the linear system
bool SchurSolver::BackSolve(double* x, double* b)
{
	// get the partition sizes
	int n0 = m_pK->PartitionEquations(0);
	int n1 = m_pK->PartitionEquations(1);

	// Get the blocks
	BlockMatrix::BLOCK& A = m_pK->Block(0, 0);
	BlockMatrix::BLOCK& B = m_pK->Block(0, 1);
	BlockMatrix::BLOCK& C = m_pK->Block(1, 0);
	BlockMatrix::BLOCK& D = m_pK->Block(1, 1);

	// split right hand side in two
	vector<double> F(n0), G(n1);
	for (int i = 0; i<n0; ++i) F[i] = m_Wu[i]*b[i];
	for (int i = 0; i<n1; ++i) G[i] = m_Wp[i]*b[i + n0];

	// solution vectors
	vector<double> u(n0, 0.0);
	vector<double> v(n1, 0.0);

	if (m_schurBlock == 0)
	{
		// step 1: solve Ay = F
		vector<double> y(n0);
		if (m_printLevel != 0) feLog("----------------------\nstep 1:\n");
		if (m_Asolver->BackSolve(y, F) == false) return false;

		// step 2: Solve Sv = H, where H = Cy - G
		if (m_printLevel != 0) feLog("step 2:\n");
		vector<double> H(n1);
		C.vmult(y, H);
		H -= G;

		if (m_schurSolver->BackSolve(v, H) == false) return false;

		// step 3: solve Au = L , where L = F - Bv
		if (m_printLevel != 0) feLog("step 3:\n");
		vector<double> tmp(n0);
		B.vmult(v, tmp);
		vector<double> L = F - tmp;
		if (m_Asolver->BackSolve(u, L) == false) return false;
	}
	else
	{
		// step 1: solve Dy = G
		vector<double> y(n1);
		if (m_printLevel != 0) feLog("----------------------\nstep 1:\n");
		if (m_Asolver->BackSolve(y, G) == false) return false;

		// step 2: Solve Su = H, where H = By - F
		if (m_printLevel != 0) feLog("step 2:\n");
		vector<double> H(n0);
		B.vmult(y, H);
		H -= F;

		if (m_schurSolver->BackSolve(u, H) == false) return false;

		// step 3: solve Dv = L , where L = G - Cu
		if (m_printLevel != 0) feLog("step 3:\n");
		vector<double> tmp(n1);
		C.vmult(u, tmp);
		vector<double> L = G - tmp;
		if (m_Asolver->BackSolve(v, L) == false) return false;
	}

	// put it back together
	for (int i = 0; i<n0; ++i) x[i     ] = m_Wu[i]*u[i];
	for (int i = 0; i<n1; ++i) x[i + n0] = m_Wp[i]*v[i];

	return true;
}

//-----------------------------------------------------------------------------
//! Clean up
void SchurSolver::Destroy()
{
	if (m_SchurAsolver != m_Asolver) m_SchurAsolver->Destroy();
	if (m_Asolver) m_Asolver->Destroy();
	if (m_schurSolver) m_schurSolver->Destroy();

	m_Asolver = nullptr;
	m_SchurAsolver = nullptr;
	m_schurSolver = nullptr;
}
