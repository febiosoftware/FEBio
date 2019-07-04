/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "NumCore.h"
#include "SkylineSolver.h"
#include "LUSolver.h"
#include "PardisoSolver.h"
#include "RCICGSolver.h"
#include "FGMRESSolver.h"
#include "FGMRES_ILU0_Solver.h"
#include "FGMRES_ILUT_Solver.h"
#include "BIPNSolver.h"
#include "HypreGMRESsolver.h"
#include "SchurSolver.h"
#include "LUPreconditioner.h"
#include "IncompleteCholesky.h"
#include "FGMRES_Schur_Solver.h"
#include "MixedLinearSolver.h"
#include "ScaledFGMRESSolver.h"
#include "FGMRES_Jacobi_Solver.h"
#include "BoomerAMGSolver.h"
#include "BlockSolver.h"
#include "BlockPreconditioner.h"
#include "FGMRES_AMG_Solver.h"
#include <FECore/fecore_enum.h>
#include <FECore/FECoreFactory.h>
#include <FECore/FECoreKernel.h>

namespace NumCore {

//=================================================================================================
class LinearSolverFactory : public FECoreFactory
{
public:
	LinearSolverFactory(const char* sztype) : FECoreFactory(FELINEARSOLVER_ID, 0, sztype) 
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterFactory(this);
	}

	LinearSolverFactory(const char* szclassName, const char* sztype) : FECoreFactory(FELINEARSOLVER_ID, szclassName, sztype)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterFactory(this);
	}
};

//=================================================================================================
class PardisoSolverFactory : public LinearSolverFactory
{
public:
	PardisoSolverFactory() : LinearSolverFactory("PardisoSolver", "pardiso")
	{
		m_estcond = false;
		m_iparm3 = false;
	}

	void* Create(FEModel* fem) const override
	{
		PardisoSolver* ls = new PardisoSolver(fem);
		ls->PrintConditionNumber(m_estcond);
		ls->UseIterativeFactorization(m_iparm3);
		return ls;
	}

private:
	bool	m_estcond;		// estimate condition number
	bool	m_iparm3;		// use iterative factorization

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(PardisoSolverFactory, FECoreFactory)
	ADD_PARAMETER(m_estcond, "print_condition_number");
	ADD_PARAMETER(m_iparm3 , "precondition");
END_FECORE_CLASS();

//=================================================================================================
class RCICGSolverFactory : public LinearSolverFactory
{
public:
	RCICGSolverFactory() : LinearSolverFactory("rcicg")
	{
		m_maxiter = 0;
		m_tol = 1e-5;
		m_print_level = 0;
	}

	void* Create(FEModel* fem) const override
	{
		RCICGSolver* ls = new RCICGSolver(fem);
		ls->SetMaxIterations(m_maxiter);
		ls->SetTolerance(m_tol);
		ls->SetPrintLevel(m_print_level);
		return ls;
	}

private:
	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_print_level;	// output level

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(RCICGSolverFactory, FECoreFactory)
	ADD_PARAMETER(m_maxiter, "maxiter");
	ADD_PARAMETER(m_tol, "tol");
	ADD_PARAMETER(m_print_level, "print_level");
END_FECORE_CLASS();

//=================================================================================================
class RCICG_ICHOL_SolverFactory : public LinearSolverFactory
{
public:
	RCICG_ICHOL_SolverFactory() : LinearSolverFactory("rcicg_ichol")
	{
		m_maxiter = 0;
		m_tol = 1e-5;
		m_print_level = 0;
	}

	void* Create(FEModel* fem) const override
	{
		RCICG_ICHOL_Solver* ls = new RCICG_ICHOL_Solver(fem);
		ls->SetMaxIterations(m_maxiter);
		ls->SetTolerance(m_tol);
		ls->SetPrintLevel(m_print_level);
		return ls;
	}

private:
	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_print_level;	// output level

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(RCICG_ICHOL_SolverFactory, FECoreFactory)
	ADD_PARAMETER(m_maxiter, "maxiter");
	ADD_PARAMETER(m_tol, "tol");
	ADD_PARAMETER(m_print_level, "print_level");
END_FECORE_CLASS();


//=================================================================================================
class FGMRES_ILUT_Factory : public LinearSolverFactory
{
public:
	FGMRES_ILUT_Factory() : LinearSolverFactory("fgmres_ilut")
	{
		m_fillTol = 1e-16;
		m_maxfill = 1;
		m_maxiter = 0; // use default min(N, 150)
		m_nrestart = 0;
		m_print_level = 0;
		m_doResidualTest = true;
		m_tol = 0.0;

		m_checkZeroDiagonal = true;
		m_zeroThreshold = 1e-16;
		m_zeroReplace = 1e-10;

		m_print_cn = false;
	}

	void* Create(FEModel* fem) const override
	{ 
		FGMRES_ILUT_Solver* ls = new FGMRES_ILUT_Solver(fem);
		ls->SetMaxFill(m_maxfill);
		ls->SetFillTolerance(m_fillTol);
		ls->SetMaxIterations(m_maxiter);
		ls->SetNonRestartedIterations(m_nrestart);
		ls->SetPrintLevel(m_print_level);
		ls->DoResidualStoppingTest(m_doResidualTest);
		ls->SetRelativeResidualTolerance(m_tol);
		ls->PrintConditionNumber(m_print_cn);

		ls->DoZeroDiagonalCheck(m_checkZeroDiagonal);
		ls->SetZeroDiagonalTolerance(m_zeroThreshold);
		ls->SetZeroDiagonalReplacement(m_zeroReplace);
		return ls;
	}

private:
	int		m_maxfill;		// max fill in values (I think this is in terms of bandwidth, not actual values)
	double	m_fillTol;		// tolerance for fill in criterion
	int		m_maxiter;			// max number of iterations
	int		m_nrestart;			// nr of non-restarted iterations
	int		m_print_level;		// print level
	bool	m_doResidualTest;	// residual stopping tets flag
	double	m_tol;				// residual convergence tolerance
	bool	m_print_cn;			// print the condition number

	// pre-conditioner parameters
	bool	m_checkZeroDiagonal;	// check for zero diagonals
	double	m_zeroThreshold;		// threshold for zero diagonal check
	double	m_zeroReplace;			// replacement value for zero diagonal

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FGMRES_ILUT_Factory, LinearSolverFactory)
	ADD_PARAMETER(m_maxfill, "maxfill");
	ADD_PARAMETER(m_fillTol, "filltol");
	ADD_PARAMETER(m_maxiter       , "maxiter");
	ADD_PARAMETER(m_nrestart      , "maxrestart");
	ADD_PARAMETER(m_print_level   , "print_level");
	ADD_PARAMETER(m_doResidualTest, "check_residual");
	ADD_PARAMETER(m_tol           , "tol");
	ADD_PARAMETER(m_checkZeroDiagonal, "replace_zero_diagonal");
	ADD_PARAMETER(m_zeroThreshold    , "zero_threshold");
	ADD_PARAMETER(m_zeroReplace      , "zero_replace");
	ADD_PARAMETER(m_print_cn         , "print_condition_number");
END_FECORE_CLASS();

//=================================================================================
class FGMRES_ILU0_Factory : public LinearSolverFactory
{
public:
	FGMRES_ILU0_Factory() : LinearSolverFactory("fgmres_ilu0")
	{
		m_maxiter = 0; // use default min(N, 150)
		m_nrestart = 0;
		m_print_level = 0;
		m_doResidualTest = true;
		m_tol = 0;
		m_abstol = 0;
		m_print_cn = false;

		m_checkZeroDiagonal = true;
		m_zeroThreshold = 1e-16;
		m_zeroReplace = 1e-10;
	}
	void* Create(FEModel* fem) const override
	{
		FGMRES_ILU0_Solver* ls = new FGMRES_ILU0_Solver(fem);
		ls->SetMaxIterations(m_maxiter);
		ls->SetNonRestartedIterations(m_nrestart);
		ls->SetPrintLevel(m_print_level);
		ls->DoResidualStoppingTest(m_doResidualTest);
		ls->SetRelativeResidualTolerance(m_tol);
		ls->SetAbsoluteResidualTolerance(m_abstol);
		ls->PrintConditionNumber(m_print_cn);
		ls->DoZeroDiagonalCheck(m_checkZeroDiagonal);
		ls->SetZeroDiagonalTolerance(m_zeroThreshold);
		ls->SetZeroDiagonalReplacement(m_zeroReplace);
		return ls;
	}

private:
	int		m_maxiter;			// max number of iterations
	int		m_nrestart;			// nr of non-restarted iterations
	int		m_print_level;		// print level
	bool	m_doResidualTest;	// residual stopping tets flag
	double	m_tol;				// residual convergence tolerance
	double	m_abstol;			// absolute convergence tolerance
	bool	m_print_cn;			// calculate and print condition number

	// pre-conditioner parameters
	bool	m_checkZeroDiagonal;	// check for zero diagonals
	double	m_zeroThreshold;		// threshold for zero diagonal check
	double	m_zeroReplace;			// replacement value for zero diagonal

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FGMRES_ILU0_Factory, LinearSolverFactory)
	ADD_PARAMETER(m_maxiter       , "maxiter");
	ADD_PARAMETER(m_nrestart      , "maxrestart");
	ADD_PARAMETER(m_print_level   , "print_level");
	ADD_PARAMETER(m_doResidualTest, "check_residual");
	ADD_PARAMETER(m_tol           , "tol");
	ADD_PARAMETER(m_abstol        , "abstol");
	ADD_PARAMETER(m_checkZeroDiagonal, "replace_zero_diagonal");
	ADD_PARAMETER(m_zeroThreshold    , "zero_threshold");
	ADD_PARAMETER(m_zeroReplace      , "zero_replace");
	ADD_PARAMETER(m_print_cn         , "print_condition_number");
END_FECORE_CLASS();

//=================================================================================
class FGMRES_Jacobi_Factory : public LinearSolverFactory
{
public:
	FGMRES_Jacobi_Factory() : LinearSolverFactory("fgmres_jacobi")
	{
		m_maxiter = 0; // use default min(N, 150)
		m_nrestart = 0;
		m_print_level = 0;
		m_doResidualTest = true;
		m_tol = 0;
		m_abstol = 0;
		m_print_cn = false;
	}

	void* Create(FEModel* fem) const override
	{
		FGMRES_Jacobi_Solver* ls = new FGMRES_Jacobi_Solver(fem);
		ls->SetMaxIterations(m_maxiter);
		ls->SetNonRestartedIterations(m_nrestart);
		ls->SetPrintLevel(m_print_level);
		ls->DoResidualStoppingTest(m_doResidualTest);
		ls->SetRelativeResidualTolerance(m_tol);
		ls->SetAbsoluteResidualTolerance(m_abstol);
		ls->PrintConditionNumber(m_print_cn);
		return ls;
	}

private:
	int		m_maxiter;			// max number of iterations
	int		m_nrestart;			// nr of non-restarted iterations
	int		m_print_level;		// print level
	bool	m_doResidualTest;	// residual stopping tets flag
	double	m_tol;				// residual convergence tolerance
	double	m_abstol;			// absolute convergence tolerance
	bool	m_print_cn;			// calculate and print condition number

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FGMRES_Jacobi_Factory, LinearSolverFactory)
	ADD_PARAMETER(m_maxiter       , "maxiter");
	ADD_PARAMETER(m_nrestart      , "maxrestart");
	ADD_PARAMETER(m_print_level   , "print_level");
	ADD_PARAMETER(m_doResidualTest, "check_residual");
	ADD_PARAMETER(m_tol           , "tol");
	ADD_PARAMETER(m_abstol        , "abstol");
	ADD_PARAMETER(m_print_cn         , "print_condition_number");
END_FECORE_CLASS();


//=================================================================================
class ScaledFGMRES_Factory : public LinearSolverFactory
{
public:
	ScaledFGMRES_Factory() : LinearSolverFactory("scaled_fgmres")
	{
		m_maxiter = 0; // use default min(N, 150)
		m_nrestart = 0;
		m_print_level = 0;
		m_doResidualTest = true;
		m_tol = 0;
		m_abstol = 0;

		m_checkZeroDiagonal = true;
		m_zeroThreshold = 1e-16;
		m_zeroReplace = 1e-10;

		m_k = 1.0;
	}
	void* Create(FEModel* fem) const override
	{
		ScaledFGMRESSolver* ls = new ScaledFGMRESSolver(fem);
		ls->SetMaxIterations(m_maxiter);
		ls->SetNonRestartedIterations(m_nrestart);
		ls->SetPrintLevel(m_print_level);
		ls->DoResidualStoppingTest(m_doResidualTest);
		ls->SetRelativeResidualTolerance(m_tol);
		ls->SetAbsoluteResidualTolerance(m_abstol);
		ls->SetScaleFactor(m_k);

//		ls->DoZeroDiagonalCheck(m_checkZeroDiagonal);
//		ls->SetZeroDiagonalTolerance(m_zeroThreshold);
//		ls->SetZeroDiagonalReplacement(m_zeroReplace);
		return ls;
	}

private:
	int		m_maxiter;			// max number of iterations
	int		m_nrestart;			// nr of non-restarted iterations
	int		m_print_level;		// print level
	bool	m_doResidualTest;	// residual stopping tets flag
	double	m_tol;				// residual convergence tolerance
	double	m_abstol;			// absolute convergence tolerance

	double	m_k;	// scale factor

	// pre-conditioner parameters
	bool	m_checkZeroDiagonal;	// check for zero diagonals
	double	m_zeroThreshold;		// threshold for zero diagonal check
	double	m_zeroReplace;			// replacement value for zero diagonal

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(ScaledFGMRES_Factory, LinearSolverFactory)
	ADD_PARAMETER(m_maxiter       , "maxiter");
	ADD_PARAMETER(m_nrestart      , "maxrestart");
	ADD_PARAMETER(m_print_level   , "print_level");
	ADD_PARAMETER(m_doResidualTest, "check_residual");
	ADD_PARAMETER(m_tol           , "tol");
	ADD_PARAMETER(m_abstol        , "abstol");
	ADD_PARAMETER(m_checkZeroDiagonal, "replace_zero_diagonal");
	ADD_PARAMETER(m_zeroThreshold    , "zero_threshold");
	ADD_PARAMETER(m_zeroReplace      , "zero_replace");
	ADD_PARAMETER(m_k               , "k");
END_FECORE_CLASS();


//=================================================================================
class FGMRESSolverFactory : public LinearSolverFactory
{
public:
	FGMRESSolverFactory() : LinearSolverFactory("fgmres")
	{
		m_maxiter = 0; // use default min(N, 150)
		m_nrestart = 0; // use maxiter
		m_print_level = 0;
		m_doResidualTest = true;
		m_tol = 0;
	}
	void* Create(FEModel* fem) const override
	{
		FGMRESSolver* ls = new FGMRESSolver(fem);
		ls->SetMaxIterations(m_maxiter);
		ls->SetNonRestartedIterations(m_nrestart);
		ls->SetPrintLevel(m_print_level);
		ls->DoResidualStoppingTest(m_doResidualTest);
		ls->SetRelativeResidualTolerance(m_tol);
		return ls;
	}

private:
	int		m_maxiter;			// max number of iterations
	int		m_nrestart;			// nr of non-restarted iterations
	int		m_print_level;		// print level
	bool	m_doResidualTest;	// residual stopping tets flag
	double	m_tol;				// residual convergence tolerance

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FGMRESSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_maxiter       , "maxiter");
	ADD_PARAMETER(m_print_level   , "print_level");
	ADD_PARAMETER(m_doResidualTest, "check_residual");
	ADD_PARAMETER(m_nrestart      , "maxrestart");
	ADD_PARAMETER(m_tol           , "tol");
END_FECORE_CLASS();

//=======================================================================================

class BIPNSolverFactory : public LinearSolverFactory
{
public:
	BIPNSolverFactory() : LinearSolverFactory("bipn")
	{
		m_maxiter = 0;
		m_tol = 1e-5;
		m_print_level = 0;
		m_use_cg = true;
		m_do_jacobi = true;

		m_cg_max = 0;
		m_cg_tol = 0;
		m_cg_res = true;

		m_gmres_max = 0;
		m_gmres_tol = 0;
		m_gmres_res = true;
		m_gmres_ilu0 = false;

		m_pc_schur = 0;
	}

	void* Create(FEModel* fem) const override
	{
		BIPNSolver* ls = new BIPNSolver(fem);
		ls->SetMaxIterations(m_maxiter);
		ls->SetTolerance(m_tol);
		ls->SetPrintLevel(m_print_level);
		ls->UseConjugateGradient(m_use_cg);
		ls->DoJacobiPreconditioner(m_do_jacobi);
		ls->SetSchurPreconditioner(m_pc_schur);

		ls->SetCGParameters(m_cg_max, m_cg_tol, m_cg_res);
		ls->SetGMRESParameters(m_gmres_max, m_gmres_tol, m_gmres_res, m_gmres_ilu0);

		return ls;
	}

private:
	// BIPN parameters
	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_print_level;	// output level
	bool	m_use_cg;		// use CG for step 2 (or GMRES otherwise)
	bool	m_do_jacobi;	// Do Jacobi preconditioning
	int		m_pc_schur;		// precondition schur system

	// CG parameters
	int		m_cg_max;
	double	m_cg_tol;
	bool	m_cg_res;

	// GMRES parameters
	int		m_gmres_max;
	double	m_gmres_tol;
	bool	m_gmres_res;
	bool	m_gmres_ilu0;

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(BIPNSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_maxiter    , "maxiter"    );
	ADD_PARAMETER(m_tol        , "tol"        );
	ADD_PARAMETER(m_print_level, "print_level");
	ADD_PARAMETER(m_use_cg     , "use_cg");
	ADD_PARAMETER(m_cg_max     , "cg_maxiter" );
	ADD_PARAMETER(m_cg_tol     , "cg_tol"     );
	ADD_PARAMETER(m_cg_res     , "cg_check_residual");
	ADD_PARAMETER(m_gmres_max  , "gmres_maxiter");
	ADD_PARAMETER(m_gmres_tol  , "gmres_tol" );
	ADD_PARAMETER(m_gmres_res  , "gmres_check_residual");
	ADD_PARAMETER(m_gmres_ilu0 , "gmres_precondition");
	ADD_PARAMETER(m_do_jacobi  , "do_jacobi");
	ADD_PARAMETER(m_pc_schur   , "precondition_schur");
END_FECORE_CLASS();

//=======================================================================================

class HYPRE_FGMRES_SolverFactory : public LinearSolverFactory
{
public:
	HYPRE_FGMRES_SolverFactory() : LinearSolverFactory("hypre_gmres")
	{
		m_maxiter = 1000;
		m_tol = 1e-7;
		m_print_level = 0;
	}

	void* Create(FEModel* fem) const override
	{
		HypreGMRESsolver* ls = new HypreGMRESsolver(fem);
		ls->SetPrintLevel(m_print_level);
		ls->SetMaxIterations(m_maxiter);
		ls->SetConvergencTolerance(m_tol);
		return ls;
	}

private:
	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_print_level;	// output level

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(HYPRE_FGMRES_SolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_print_level, "print_level");
	ADD_PARAMETER(m_maxiter    , "maxiter"    );
	ADD_PARAMETER(m_tol        , "tol"        );
END_FECORE_CLASS();


//=============================================================================
class SchurLinearSolverFactory : public LinearSolverFactory
{
public:
	SchurLinearSolverFactory() : LinearSolverFactory("schur")
	{
		m_maxiter = 0;
		m_reltol = 1e-7;
		m_abstol = 0.0;
		m_print_level = 0;
		m_nAsolver = 0;
		m_nSchurSolver = 0;
		m_nSchurASolver = -1;
		m_nSchurPC = 0;
		m_doJacobi = false;
		m_schurBlock = 0;
	}

	void* Create(FEModel* fem) const override
	{
		SchurSolver* ls = new SchurSolver(fem);
		ls->SetPrintLevel(m_print_level);
		ls->SetMaxIterations(m_maxiter);
		ls->SetRelativeResidualTolerance(m_reltol);
		ls->SetAbsoluteResidualTolerance(m_abstol);
		ls->SetLinearSolver(m_nAsolver);
		ls->SetSchurSolver(m_nSchurSolver);
		ls->SetSchurPreconditioner(m_nSchurPC);
		ls->DoJacobiPreconditioning(m_doJacobi);
		ls->SetSchurBlock(m_schurBlock);
		ls->SetSchurASolver(m_nSchurASolver == -1 ? m_nAsolver : m_nSchurASolver);
		return ls;
	}

private:
	int		m_maxiter;		// max nr of iterations
	double	m_reltol;		// residual relative tolerance
	double	m_abstol;		// residual absolute tolerance
	int		m_print_level;	// output level
	int		m_nAsolver;		// A block solver
	int		m_nSchurSolver;	// Schur complement solver
	int		m_nSchurASolver;// Schur A-block solver
	int		m_nSchurPC;		// Schur complement preconditioner
	int		m_schurBlock;	// which block to use for Schur solver
	bool	m_doJacobi;		// Do Jacobi preconditioning

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(SchurLinearSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_print_level , "print_level");
	ADD_PARAMETER(m_maxiter     , "maxiter");
	ADD_PARAMETER(m_reltol      , "tol");
	ADD_PARAMETER(m_abstol      , "abstol");
	ADD_PARAMETER(m_nAsolver    , "linear_solver");
	ADD_PARAMETER(m_nSchurSolver, "schur_solver");
	ADD_PARAMETER(m_nSchurASolver, "schur_Asolver");
	ADD_PARAMETER(m_nSchurPC    , "schur_pc");
	ADD_PARAMETER(m_doJacobi    , "do_jacobi");
	ADD_PARAMETER(m_schurBlock  , "schur_block");
END_FECORE_CLASS();

//=============================================================================
class FGMRESSchurLinearSolverFactory : public LinearSolverFactory
{
public:
	FGMRESSchurLinearSolverFactory() : LinearSolverFactory("fgmres_schur")
	{
		m_maxiter = 0;
		m_print_level = 0;
		m_tol = 1e-6;
		m_bzeroDBlock = false;

		m_do_jacobi = false;

		m_solver = 0;
		m_schurSolver = 2;
		m_schurPC = 1;
		m_schurBlock = 0;
		m_schurASolver = -1;

		m_schur_maxiter = 0;
		m_schur_tol = 1e-8;
	}

	void* Create(FEModel* fem) const override
	{
		FGMRES_Schur_Solver* ls = new FGMRES_Schur_Solver(fem);
		ls->SetPrintLevel(m_print_level);
		ls->SetMaxIterations(m_maxiter);
		ls->SetRelativeResidualTolerance(m_tol);
		ls->ZeroDBlock(m_bzeroDBlock);

		SchurPreconditioner* pc = dynamic_cast<SchurPreconditioner*>(ls->GetPreconditioner());
		pc->SetLinearSolver(m_solver);
		pc->SetSchurSolver(m_schurSolver);
		pc->SetSchurPreconditioner(m_schurPC);
		pc->SetMaxIterations(m_schur_maxiter);
		pc->SetTolerance(m_schur_tol);
		pc->SetSchurBlock(m_schurBlock);
		pc->SetSchurASolver(m_schurASolver == -1 ? m_solver : m_schurASolver);
		pc->DoJacobiPreconditioning(m_do_jacobi);
		return ls;
	}

private:
	int		m_maxiter;		// max nr of iterations
	int		m_print_level;	// output level
	double	m_tol;
	bool	m_bzeroDBlock;

	bool	m_do_jacobi;

	int		m_schur_maxiter;
	double	m_schur_tol;

	int m_solver;
	int	m_schurSolver;
	int m_schurPC;
	int m_schurBlock;
	int m_schurASolver;

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FGMRESSchurLinearSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_print_level  , "print_level");
	ADD_PARAMETER(m_maxiter      , "maxiter");
	ADD_PARAMETER(m_tol          , "tol");
	ADD_PARAMETER(m_bzeroDBlock  , "zero_D_block");
	ADD_PARAMETER(m_solver       , "linear_solver");
	ADD_PARAMETER(m_schurSolver  , "schur_solver");
	ADD_PARAMETER(m_schurPC      , "schur_pc");
	ADD_PARAMETER(m_schur_maxiter, "schur_maxiter");
	ADD_PARAMETER(m_schur_tol    , "schur_tol");
	ADD_PARAMETER(m_schurASolver , "schur_Asolver");
	ADD_PARAMETER(m_schurBlock   , "schur_block");
	ADD_PARAMETER(m_do_jacobi    , "do_jacobi");
END_FECORE_CLASS();

//=================================================================================
class MixedSolverFactory : public LinearSolverFactory
{
public:
	MixedSolverFactory() : LinearSolverFactory("mixed")
	{
		m_maxiter = 0; // use default min(N, 150)
		m_print_level = 0;
		m_relTol = 0;
		m_absTol = 0;
	}

	void* Create(FEModel* fem) const override
	{
		MixedLinearSolver* ls = new MixedLinearSolver(fem);
		ls->SetMaxIterations(m_maxiter);
		ls->SetPrintLevel(m_print_level);
		ls->SetRelativeConvergence(m_relTol);
		ls->SetAbsoluteConvergence(m_absTol);
		return ls;
	}

private:
	int		m_maxiter;			// max number of iterations
	int		m_print_level;		// print level
	double	m_relTol;			// residual convergence tolerance
	double	m_absTol;			// absolute convergence tolerance

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(MixedSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_maxiter       , "maxiter");
	ADD_PARAMETER(m_print_level   , "print_level");
	ADD_PARAMETER(m_relTol        , "tol");
	ADD_PARAMETER(m_absTol        , "abstol");
END_FECORE_CLASS();

//=================================================================================
class BoomerAMGSolverFactory : public LinearSolverFactory
{
public:
	BoomerAMGSolverFactory() : LinearSolverFactory("BoomerAMGSolver", "boomeramg")
	{
		m_maxiter = 0; // use default min(N, 150)
		m_print_level = 0;
		m_relTol = 0;
		m_maxLevels = 25;
		m_coarsenType = -1;
		m_use_num_funcs = false;
        m_relaxType = 3;    /* hybrid Gauss-Seidel or SOR, forward solve */
        m_interpType = 6;   /* extended+i interpolation */
        m_strong_threshold = 0.5;
        m_PMaxElmts = 4;
        m_NumSweeps = 1;
        m_AggNumLevels = 0;
		m_nodal = 0;
	}

	void* Create(FEModel* fem) const override
	{
		BoomerAMGSolver* ls = new BoomerAMGSolver(fem);
		ls->SetMaxIterations(m_maxiter);
		ls->SetPrintLevel(m_print_level);
		ls->SetConvergenceTolerance(m_relTol);
		ls->SetMaxLevels(m_maxLevels);
		ls->SetCoarsenType(m_coarsenType);
		ls->SetUseNumFunctions(m_use_num_funcs);
        ls->SetRelaxType(m_relaxType);
        ls->SetInterpType(m_interpType);
        ls->SetStrongThreshold(m_strong_threshold);
        ls->SetPMaxElmts(m_PMaxElmts);
        ls->SetNumSweeps(m_NumSweeps);
        ls->SetAggNumLevels(m_AggNumLevels);
		ls->SetNodal(m_nodal);
		return ls;
	}

private:
	int		m_maxiter;			// max number of iterations
	int		m_maxLevels;		// set max of multigrid levels
	int		m_print_level;		// print level
	double	m_relTol;			// residual convergence tolerance
	int		m_coarsenType;		// set coarsening type
	bool	m_use_num_funcs;	// use the number of functions feature
    int     m_relaxType;
    int     m_interpType;
    int     m_PMaxElmts;
    int     m_NumSweeps;
    int     m_AggNumLevels;
    double  m_strong_threshold;
	int		m_nodal;

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(BoomerAMGSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_maxiter         , "max_iter");
	ADD_PARAMETER(m_print_level     , "print_level");
	ADD_PARAMETER(m_relTol          , "tol");
	ADD_PARAMETER(m_maxLevels       , "max_levels");
	ADD_PARAMETER(m_coarsenType     , "coarsen_type");
	ADD_PARAMETER(m_use_num_funcs	, "use_num_funcs");
    ADD_PARAMETER(m_relaxType       , "relax_type");
    ADD_PARAMETER(m_interpType      , "interp_type");
    ADD_PARAMETER(m_strong_threshold, "strong_threshold");
    ADD_PARAMETER(m_PMaxElmts       , "p_max_elmts");
    ADD_PARAMETER(m_NumSweeps       , "num_sweeps");
    ADD_PARAMETER(m_AggNumLevels    , "agg_num_levels");
	ADD_PARAMETER(m_nodal           , "nodal");
END_FECORE_CLASS();


//=================================================================================
class BlockJacobiSolverFactory : public LinearSolverFactory
{
public:
	BlockJacobiSolverFactory() : LinearSolverFactory("block-jacobi")
	{
		m_relTol = 1e-8;
		m_maxIter = 150;
		m_printLevel = 0;
		m_method = BlockIterativeSolver::JACOBI;
		m_failMaxIter = true;
		m_zeroInitGuess = true;
	}

	void* Create(FEModel* fem) const override
	{
		BlockIterativeSolver* ls = new BlockIterativeSolver(fem);

		ls->SetRelativeTolerance(m_relTol);
		ls->SetMaxIterations(m_maxIter);
		ls->SetPrintLevel(m_printLevel);
		ls->SetFailMaxIters(m_failMaxIter);
		ls->SetSolutionMethod(m_method);
		ls->SetZeroInitialGuess(m_zeroInitGuess);

		return ls;
	}

private:
	double	m_relTol;		// relative tolerance
	int		m_maxIter;		// max number of iterations
	int		m_printLevel;	// output level
	bool	m_failMaxIter;		//!< fail on max iterations reached
	int		m_method;			//!< set solution method
	bool	m_zeroInitGuess;	//!< always use zero as the initial guess

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(BlockJacobiSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_maxIter   , "max_iter");
	ADD_PARAMETER(m_printLevel, "print_level");
	ADD_PARAMETER(m_relTol    , "tol");
	ADD_PARAMETER(m_failMaxIter, "fail_max_iter");
	ADD_PARAMETER(m_method     , "solution_method");
	ADD_PARAMETER(m_zeroInitGuess, "zero_initial_guess");
END_FECORE_CLASS()


//=================================================================================
class FGMRES_Jacobi_Block_Solver_Factory : public LinearSolverFactory
{
public:
	FGMRES_Jacobi_Block_Solver_Factory() : LinearSolverFactory("fgmres_block-jacobi")
	{
		m_relTol = 1e-8;
		m_maxIter = 150;
		m_printLevel = 0;
		m_method = 0;
		m_blockSolver = 0; // pardiso
	}

	void* Create(FEModel* fem) const override
	{
		FGMRES_Jacobi_Block_Solver* ls = new FGMRES_Jacobi_Block_Solver(fem);

		ls->SetRelativeResidualTolerance(m_relTol);
		ls->SetMaxIterations(m_maxIter);
		ls->SetPrintLevel(m_printLevel);
		ls->SetSolutionMethod(m_method);
		ls->SetBlockSolver(m_blockSolver);

		return ls;
	}

private:
	double	m_relTol;		// relative tolerance
	int		m_maxIter;		// max number of iterations
	int		m_printLevel;	// output level
	int		m_method;		//!< set solution method
	int		m_blockSolver;	//!< block solver used by preconditioner

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FGMRES_Jacobi_Block_Solver_Factory, LinearSolverFactory)
	ADD_PARAMETER(m_maxIter   , "max_iter");
	ADD_PARAMETER(m_printLevel, "print_level");
	ADD_PARAMETER(m_relTol    , "tol");
	ADD_PARAMETER(m_method    , "solution_method");
	ADD_PARAMETER(m_blockSolver, "block_solver", 0, "@factory_list:31");
END_FECORE_CLASS()

//=================================================================================
class FGMRES_AMG_Factory : public LinearSolverFactory
{
public:
	FGMRES_AMG_Factory() : LinearSolverFactory("fgmres_amg")
	{
		m_relTol = 1e-8;
		m_maxIter = 150;
		m_printLevel = 0;
	}

	void* Create(FEModel* fem) const override
	{
		FGMRES_AMG_Solver* ls = new FGMRES_AMG_Solver(fem);

		ls->SetRelativeResidualTolerance(m_relTol);
		ls->SetMaxIterations(m_maxIter);
		ls->SetPrintLevel(m_printLevel);

		return ls;
	}

private:
	double	m_relTol;		// relative tolerance
	int		m_maxIter;		// max number of iterations
	int		m_printLevel;	// output level

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FGMRES_AMG_Factory, LinearSolverFactory)
	ADD_PARAMETER(m_maxIter   , "max_iter");
	ADD_PARAMETER(m_printLevel, "print_level");
	ADD_PARAMETER(m_relTol    , "tol");
END_FECORE_CLASS()

} // namespace NumCore

//=============================================================================
// Call this to initialize the NumCore module
void NumCore::InitModule()
{
	// register linear solvers
	REGISTER_FECORE_FACTORY(PardisoSolverFactory          );
	REGISTER_FECORE_FACTORY(SchurLinearSolverFactory      );
	REGISTER_FECORE_FACTORY(HYPRE_FGMRES_SolverFactory    );
	REGISTER_FECORE_FACTORY(BIPNSolverFactory             );
	REGISTER_FECORE_FACTORY(FGMRESSolverFactory           );
	REGISTER_FECORE_FACTORY(FGMRES_ILU0_Factory           );
	REGISTER_FECORE_FACTORY(FGMRES_ILUT_Factory           );
	REGISTER_FECORE_FACTORY(RCICGSolverFactory            );
	REGISTER_FECORE_FACTORY(RCICG_ICHOL_SolverFactory     );
	REGISTER_FECORE_FACTORY(FGMRESSchurLinearSolverFactory);
	REGISTER_FECORE_FACTORY(MixedSolverFactory            );
	REGISTER_FECORE_FACTORY(ScaledFGMRES_Factory          );
	REGISTER_FECORE_FACTORY(FGMRES_Jacobi_Factory         );
	REGISTER_FECORE_FACTORY(BoomerAMGSolverFactory        );
	REGISTER_FECORE_FACTORY(BlockJacobiSolverFactory      );
	REGISTER_FECORE_FACTORY(FGMRES_Jacobi_Block_Solver_Factory);
	REGISTER_FECORE_FACTORY(FGMRES_AMG_Factory            );

	REGISTER_FECORE_CLASS(SkylineSolver    , "skyline"  );
	REGISTER_FECORE_CLASS(LUSolver         , "LU"       );

	// register preconditioners
	REGISTER_FECORE_CLASS(ILU0_Preconditioner, "ilu0");
	REGISTER_FECORE_CLASS(ILUT_Preconditioner, "ilut");
	REGISTER_FECORE_CLASS(LUPreconditioner   , "pardiso");
	REGISTER_FECORE_CLASS(IncompleteCholesky , "ichol");

	// set default linear solver
	// (Set this before the configuration is read in because
	//  the configuration can change the default linear solver.)
	FECoreKernel& fecore = FECoreKernel::GetInstance();
#ifdef PARDISO
	fecore.SetDefaultSolver("pardiso");
#else
	fecore.SetDefaultSolver("skyline");
#endif
}
