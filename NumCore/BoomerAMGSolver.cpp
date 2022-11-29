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
#include "BoomerAMGSolver.h"
#include <FECore/CompactUnSymmMatrix.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FESolver.h>
#include <FECore/FEAnalysis.h>
#ifdef HYPRE
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_utilities.h>
#include <_hypre_IJ_mv.h>

class BoomerAMGSolver::Implementation
{
public:
	FEModel*			m_fem;

	CRSSparseMatrix*	m_A;
	HYPRE_IJMatrix		ij_A;
	HYPRE_ParCSRMatrix	par_A;

	vector<int>			m_ind;
	HYPRE_IJVector		ij_b, ij_x;
	HYPRE_ParVector		par_b, par_x;

	vector<double>		W;		// Jacobi preconditioner
	vector<double>		rhs;	// right-hand side 

	HYPRE_Solver	m_solver;
	HYPRE_Int		m_num_iterations;
	double			m_final_res_norm;

	int			m_print_level;
	int			m_maxIter;
	int			m_maxLevels;
	double		m_tol;
	int			m_coarsenType;
    bool        m_set_num_funcs;
    int         m_relaxType;
    int         m_interpType;
    int         m_PMaxElmts;
    int         m_NumSweeps;
    int         m_AggInterpType;
    int         m_AggNumLevels;
    double      m_strong_threshold;
	int			m_nodal;
	bool		m_jacobi_pc;
	bool		m_failMaxIters;

	HYPRE_Int*	m_dofMap;

public:
	Implementation()
	{
		m_print_level = 0;
		m_maxLevels = 25;
		m_maxIter = 20;
		m_tol = 1.e-7;
		m_coarsenType = -1;	// don't set
		m_set_num_funcs = false;
        m_relaxType = 3;
        m_interpType = 6;
        m_strong_threshold = 0.5;
        m_PMaxElmts = 4;
        m_NumSweeps = 1;
        m_AggInterpType = 0;
        m_AggNumLevels = 0;
		m_nodal = 0;
		m_jacobi_pc = false;
		m_failMaxIters = true;

		m_A = nullptr;
		m_solver = nullptr;
		ij_A = nullptr;
		ij_b = nullptr;
		ij_x = nullptr;
	}

	int equations() const { return (m_A ? m_A->Rows() : 0); }

	// Allocate stiffness matrix
	void allocMatrix()
	{
		int neq = equations();

		rhs.assign(neq, 0.0);

		// Create an empty matrix object
		int ret = 0;
		ret = HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, neq - 1, 0, neq - 1, &ij_A);

		// set the matrix object type
		ret = HYPRE_IJMatrixSetObjectType(ij_A, HYPRE_PARCSR);
	}

	// destroy stiffness matrix
	void destroyMatrix()
	{
		if (ij_A) HYPRE_IJMatrixDestroy(ij_A);
		ij_A = nullptr;
	}

	// update coefficient matrix
	bool updateMatrix()
	{
		int neq = equations();

		SparseMatrix* A = m_A;
		W.resize(neq, 1.0);
		if (m_jacobi_pc)
		{
			for (int i = 0; i < neq; ++i)
			{
				double dii = fabs(A->diag(i));
				if (dii == 0.0) return false;
				W[i] = 1.0 / sqrt(dii);
			}
		}

		if (m_jacobi_pc)
			A->scale(W, W);

		// call initialize, after which we can set the matrix coefficients
		HYPRE_Int ret = HYPRE_IJMatrixInitialize(ij_A);

		// set the matrix coefficients
		double* values = m_A->Values();
		int* indices = m_A->Indices();
		int* pointers = m_A->Pointers();
		for (int i = 0; i<neq; ++i)
		{
			const int* cols = indices + pointers[i];
			int ncols = pointers[i + 1] - pointers[i];
			double* vals = values + pointers[i];
			HYPRE_Int nrows = 1;
			ret = HYPRE_IJMatrixSetValues(ij_A, nrows, (HYPRE_Int*)&ncols, (HYPRE_Int*)&i, (HYPRE_Int*)cols, vals);
		}

		// Finalize the matrix assembly
		ret = HYPRE_IJMatrixAssemble(ij_A);

		// get the matrix object for later use
		ret = HYPRE_IJMatrixGetObject(ij_A, (void**)&par_A);

		return true;
	}

	// Allocate vectors for rhs and solution
	void allocVectors()
	{
		int neq = equations();
		m_ind.resize(neq, 0);
		for (int i = 0; i<neq; ++i) m_ind[i] = i;

		// create the vector object for the rhs
		HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, neq - 1, &ij_b);
		HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);

		// create the vector object for the solution
		HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, neq - 1, &ij_x);
		HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
	}

	// destroy vectors
	void destroyVectors()
	{
		if (ij_b) HYPRE_IJVectorDestroy(ij_b); ij_b = nullptr;
		if (ij_x) HYPRE_IJVectorDestroy(ij_x); ij_x = nullptr;
	}

	// update vectors 
	void updateVectors(double* x, double* b)
	{
		// initialize vectors for changing coefficient values
		HYPRE_IJVectorInitialize(ij_b);
		HYPRE_IJVectorInitialize(ij_x);

		// scale by Jacobi PC
		int neq = equations();
		for (int i = 0; i < neq; ++i) rhs[i] = W[i] * b[i];

		// set the values
		HYPRE_IJVectorSetValues(ij_b, neq, (HYPRE_Int*)&m_ind[0], &rhs[0]);
		HYPRE_IJVectorSetValues(ij_x, neq, (HYPRE_Int*)&m_ind[0], x);

		// finialize assembly
		HYPRE_IJVectorAssemble(ij_b);
		HYPRE_IJVectorAssemble(ij_x);

		HYPRE_IJVectorGetObject(ij_b, (void**)&par_b);
		HYPRE_IJVectorGetObject(ij_x, (void**)&par_x);
	}

	bool allocSolver()
	{
		// create solver
		HYPRE_BoomerAMGCreate(&m_solver);

		int printLevel = m_print_level;
		if (printLevel > 3) printLevel = 0;

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_BoomerAMGSetPrintLevel(m_solver, printLevel);  // print solve info
//		HYPRE_BoomerAMGSetOldDefault(m_solver); /* Falgout coarsening with modified classical interpolation */
		if (m_coarsenType >= 0) HYPRE_BoomerAMGSetCoarsenType(m_solver, m_coarsenType);
		if (m_relaxType   >= 0) HYPRE_BoomerAMGSetRelaxType(m_solver, m_relaxType);
		HYPRE_BoomerAMGSetRelaxOrder(m_solver, 1);   // uses C/F relaxation
        HYPRE_BoomerAMGSetInterpType(m_solver, m_interpType);
        HYPRE_BoomerAMGSetPMaxElmts(m_solver, m_PMaxElmts);
		HYPRE_BoomerAMGSetNumSweeps(m_solver, m_NumSweeps);   /* Sweeps on each level */
		HYPRE_BoomerAMGSetMaxLevels(m_solver, m_maxLevels);  /* maximum number of levels */
        HYPRE_BoomerAMGSetStrongThreshold(m_solver, m_strong_threshold);
        HYPRE_BoomerAMGSetAggInterpType(m_solver, m_AggInterpType);
        HYPRE_BoomerAMGSetAggNumLevels(m_solver, m_AggNumLevels);
		HYPRE_BoomerAMGSetMaxIter(m_solver, m_maxIter);
		HYPRE_BoomerAMGSetTol(m_solver, m_tol);      /* conv. tolerance */

		// sets the number of funtions, which I think means the dofs per variable
		// For 3D mechanics this would be 3.
		// I think this also requires the staggered equation numbering.
		// NOTE: when set to 3, this definitely seems to help with (some) mechanics problems!
		if (m_set_num_funcs)
		{
			// we need the FESolver so we can get the dof mapping
			FESolver* fesolve = m_fem->GetCurrentStep()->GetFESolver();

			// get the dof map
			vector<int> dofMap;
			int nfunc = fesolve->GetActiveDofMap(dofMap);
			if (nfunc == -1) return false;

			// allocate dof map
			// (We need to copy it here since Hypre will deallocate it)
			int neq = (int)dofMap.size();
			m_dofMap = (HYPRE_Int*)malloc(neq * sizeof(HYPRE_Int));
			for (size_t i = 0; i < neq; ++i) m_dofMap[i] = dofMap[i];

			printf("\tNumber of functions : %d\n", nfunc);

			// assign to BoomerAMG
			HYPRE_BoomerAMGSetNumFunctions(m_solver, nfunc);

			// set the dof map
			HYPRE_BoomerAMGSetDofFunc(m_solver, m_dofMap);
		}

		// NOTE: Turning this option on seems to worsen convergence!
		if (m_nodal > 0) HYPRE_BoomerAMGSetNodal(m_solver, m_nodal);

		return true;
	}

	void destroySolver()
	{
		if (m_solver) HYPRE_BoomerAMGDestroy(m_solver);
		m_solver = nullptr;
	}

	void doSetup()
	{
		HYPRE_BoomerAMGSetup(m_solver, par_A, par_b, par_x);
	}

	bool doSolve(double* x)
	{
		HYPRE_BoomerAMGSolve(m_solver, par_A, par_b, par_x);

		/* Run info - needed logging turned on */
		HYPRE_BoomerAMGGetNumIterations(m_solver, &m_num_iterations);
		HYPRE_BoomerAMGGetFinalRelativeResidualNorm(m_solver, &m_final_res_norm);

		// see if we converged
		bool bok = false;
		if ((m_num_iterations <= m_maxIter) && (m_final_res_norm < m_tol)) bok = true;

		/* get the local solution */
		int neq = equations();
		HYPRE_IJVectorGetValues(ij_x, neq, (HYPRE_Int*)&m_ind[0], &x[0]);

		// scale by Jacobi PC
		for (int i = 0; i < neq; ++i) x[i] *= W[i];

		return bok;
	}
};

BEGIN_FECORE_CLASS(BoomerAMGSolver, LinearSolver)
	ADD_PARAMETER(imp->m_maxIter         , "max_iter");
	ADD_PARAMETER(imp->m_print_level     , "print_level");
	ADD_PARAMETER(imp->m_tol             , "tol");
	ADD_PARAMETER(imp->m_maxLevels       , "max_levels");
	ADD_PARAMETER(imp->m_coarsenType     , "coarsen_type");
	ADD_PARAMETER(imp->m_set_num_funcs	 , "use_num_funcs");
    ADD_PARAMETER(imp->m_relaxType       , "relax_type");
    ADD_PARAMETER(imp->m_interpType      , "interp_type");
    ADD_PARAMETER(imp->m_strong_threshold, "strong_threshold");
    ADD_PARAMETER(imp->m_PMaxElmts       , "p_max_elmts");
    ADD_PARAMETER(imp->m_NumSweeps       , "num_sweeps");
    ADD_PARAMETER(imp->m_AggInterpType   , "agg_interp_type");
    ADD_PARAMETER(imp->m_AggNumLevels    , "agg_num_levels");
	ADD_PARAMETER(imp->m_nodal           , "nodal");
	ADD_PARAMETER(imp->m_jacobi_pc       , "do_jacobi");
	ADD_PARAMETER(imp->m_failMaxIters    , "fail_max_iters");
END_FECORE_CLASS();

BoomerAMGSolver::BoomerAMGSolver(FEModel* fem) : LinearSolver(fem), imp(new BoomerAMGSolver::Implementation)
{
	imp->m_fem = fem;
}

BoomerAMGSolver::~BoomerAMGSolver()
{
	Destroy();
	delete imp;
}

void BoomerAMGSolver::SetPrintLevel(int printLevel)
{
	imp->m_print_level = printLevel;
}

void BoomerAMGSolver::SetMaxIterations(int maxIter)
{
	imp->m_maxIter = maxIter;
}

void BoomerAMGSolver::SetMaxLevels(int levels)
{
	imp->m_maxLevels = levels;
}

void BoomerAMGSolver::SetConvergenceTolerance(double tol)
{
	imp->m_tol = tol;
}

void BoomerAMGSolver::SetCoarsenType(int coarsenType)
{
	imp->m_coarsenType = coarsenType;
}

void BoomerAMGSolver::SetUseNumFunctions(bool b)
{
	imp->m_set_num_funcs = b;
}

void BoomerAMGSolver::SetRelaxType(int rlxtyp)
{
    imp->m_relaxType = rlxtyp;
}

void BoomerAMGSolver::SetInterpType(int inptyp)
{
    imp->m_interpType = inptyp;
}

void BoomerAMGSolver::SetStrongThreshold(double thresh)
{
    imp->m_strong_threshold = thresh;
}

void BoomerAMGSolver::SetPMaxElmts(int pmax)
{
    imp->m_PMaxElmts = pmax;
}

void BoomerAMGSolver::SetNumSweeps(int nswp)
{
    imp->m_NumSweeps = nswp;
}

void BoomerAMGSolver::SetAggInterpType(int aggit)
{
    imp->m_AggInterpType = aggit;
}

void BoomerAMGSolver::SetAggNumLevels(int anlv)
{
    imp->m_AggNumLevels = anlv;
}

void BoomerAMGSolver::SetNodal(int nodal)
{
	imp->m_nodal = nodal;
}

void BoomerAMGSolver::SetJacobiPC(bool b)
{
	imp->m_jacobi_pc = b;
}

bool BoomerAMGSolver::GetJacobiPC()
{
	return imp->m_jacobi_pc;
}

void BoomerAMGSolver::SetFailOnMaxIterations(bool b)
{
	imp->m_failMaxIters = b;
}

SparseMatrix* BoomerAMGSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	// allocate the correct matrix format depending on matrix symmetry type
	switch (ntype)
	{
	case REAL_SYMMETRIC: imp->m_A = nullptr; break;
	case REAL_UNSYMMETRIC: 
	case REAL_SYMM_STRUCTURE:
		imp->m_A = new CRSSparseMatrix(0); break;
	default:
		assert(false);
		imp->m_A = nullptr;
	}

	return imp->m_A;
}

bool BoomerAMGSolver::SetSparseMatrix(SparseMatrix* pA)
{
	if (pA == imp->m_A) return true;

	if (imp->m_A) delete imp->m_A;
	imp->m_A = dynamic_cast<CRSSparseMatrix*>(pA);
	if (imp->m_A == nullptr) return false;
	if (imp->m_A->Offset() != 0)
	{
		imp->m_A = nullptr;
		return false;
	}
	return true;
}

bool BoomerAMGSolver::PreProcess()
{
	imp->m_fem = GetFEModel();
	imp->allocMatrix();
	imp->allocVectors();

	// create solver
	return imp->allocSolver();
}

bool BoomerAMGSolver::Factor()
{
	if (imp->updateMatrix() == false) return false;

	vector<double> zero(imp->equations(), 0.0);
	imp->updateVectors(&zero[0], &zero[0]);

	imp->doSetup();

	return true;
}

bool BoomerAMGSolver::BackSolve(double* x, double* b)
{
	imp->updateVectors(x, b);

	bool bok = imp->doSolve(x);

	if (imp->m_print_level > 3)
	{
		feLog("AMG: %d iterations, %lg residual norm\n", imp->m_num_iterations, imp->m_final_res_norm);
	}

	UpdateStats(imp->m_num_iterations);

	return (bok || (imp->m_failMaxIters == false));
}

void BoomerAMGSolver::Destroy()
{
	// Destroy solver
	imp->destroySolver();

	imp->destroyMatrix();
	imp->destroyVectors();
}

#else
BEGIN_FECORE_CLASS(BoomerAMGSolver, LinearSolver)
END_FECORE_CLASS();
BoomerAMGSolver::BoomerAMGSolver(FEModel* fem) : LinearSolver(fem) {}
BoomerAMGSolver::~BoomerAMGSolver() {}
void BoomerAMGSolver::SetPrintLevel(int printLevel) {}
void BoomerAMGSolver::SetMaxIterations(int maxIter) {}
void BoomerAMGSolver::SetConvergenceTolerance(double tol) {}
void BoomerAMGSolver::SetMaxLevels(int levels) {}
void BoomerAMGSolver::SetCoarsenType(int coarsenType) {}
void BoomerAMGSolver::SetRelaxType(int rlxtyp) {}
void BoomerAMGSolver::SetInterpType(int inptyp) {}
void BoomerAMGSolver::SetStrongThreshold(double thresh) {}
void BoomerAMGSolver::SetPMaxElmts(int pmax) {}
void BoomerAMGSolver::SetNumSweeps(int nswp) {}
void BoomerAMGSolver::SetAggInterpType(int aggit) {}
void BoomerAMGSolver::SetAggNumLevels(int anlv) {}
SparseMatrix* BoomerAMGSolver::CreateSparseMatrix(Matrix_Type ntype) { return nullptr; }
bool BoomerAMGSolver::SetSparseMatrix(SparseMatrix* pA) { return false; }
bool BoomerAMGSolver::PreProcess() { return false; }
bool BoomerAMGSolver::Factor() { return false; }
bool BoomerAMGSolver::BackSolve(double* x, double* y) { return false; }
void BoomerAMGSolver::Destroy() {}
#endif
