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
#include <FECore/fecore_enum.h>
#include <FECore/FECoreFactory.h>
#include <FECore/FECoreKernel.h>

namespace NumCore {

//=================================================================================================
class LinearSolverFactory : public FECoreFactory
{
public:
	LinearSolverFactory(const char* sztype) : FECoreFactory(FELINEARSOLVER_ID, sztype) 
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterFactory(this);
	}
};

//=================================================================================================
class PardisoSolverFactory : public LinearSolverFactory
{
public:
	PardisoSolverFactory() : LinearSolverFactory("pardiso")
	{
		m_estcond = false;
	}

	void* Create(FEModel* fem) override
	{
		PardisoSolver* ls = new PardisoSolver(fem);
		ls->PrintConditionNumber(m_estcond);
		return ls;
	}

private:
	bool	m_estcond;		// estimate condition number

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(PardisoSolverFactory, FECoreFactory)
	ADD_PARAMETER(m_estcond, "print_condition_number");
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

	void* Create(FEModel* fem) override
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

	void* Create(FEModel* fem) override
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

	}
	void* Create(FEModel* fem) override
	{ 
		FGMRES_ILUT_Solver* ls = new FGMRES_ILUT_Solver(fem);
		ls->SetMaxFill(m_maxfill);
		ls->SetFillTolerance(m_fillTol);
		ls->SetMaxIterations(m_maxiter);
		ls->SetNonRestartedIterations(m_nrestart);
		ls->SetPrintLevel(m_print_level);
		ls->DoResidualStoppingTest(m_doResidualTest);
		ls->SetRelativeResidualTolerance(m_tol);

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

		m_checkZeroDiagonal = true;
		m_zeroThreshold = 1e-16;
		m_zeroReplace = 1e-10;
	}
	void* Create(FEModel* fem) override
	{
		FGMRES_ILU0_Solver* ls = new FGMRES_ILU0_Solver(fem);
		ls->SetMaxIterations(m_maxiter);
		ls->SetNonRestartedIterations(m_nrestart);
		ls->SetPrintLevel(m_print_level);
		ls->DoResidualStoppingTest(m_doResidualTest);
		ls->SetRelativeResidualTolerance(m_tol);
		ls->SetAbsoluteResidualTolerance(m_abstol);

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
	void* Create(FEModel* fem) override
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
	void* Create(FEModel* fem) override
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

		m_cg_max = 0;
		m_cg_tol = 0;
		m_cg_res = true;

		m_gmres_max = 0;
		m_gmres_tol = 0;
		m_gmres_res = true;
		m_gmres_ilu0 = false;
	}

	void* Create(FEModel* fem) override
	{
		BIPNSolver* ls = new BIPNSolver(fem);
		ls->SetMaxIterations(m_maxiter);
		ls->SetTolerance(m_tol);
		ls->SetPrintLevel(m_print_level);
		ls->UseConjugateGradient(m_use_cg);

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

	void* Create(FEModel* fem) override
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
		m_nsolver = 0;
		m_nschurSolver = 0;
		m_k = 1.0;
	}

	void* Create(FEModel* fem) override
	{
		SchurSolver* ls = new SchurSolver(fem);
		ls->SetPrintLevel(m_print_level);
		ls->SetMaxIterations(m_maxiter);
		ls->SetRelativeResidualTolerance(m_reltol);
		ls->SetAbsoluteResidualTolerance(m_abstol);
		ls->SetLinearSolver(m_nsolver);
		ls->SetSchurSolver(m_nschurSolver);
		ls->SetScaleFactor(m_k);
		return ls;
	}

private:
	int		m_maxiter;		// max nr of iterations
	double	m_reltol;		// residual relative tolerance
	double	m_abstol;		// residual absolute tolerance
	int		m_print_level;	// output level
	int		m_nsolver;
	int		m_nschurSolver;

	double	m_k;

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(SchurLinearSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_print_level , "print_level");
	ADD_PARAMETER(m_maxiter     , "maxiter");
	ADD_PARAMETER(m_reltol      , "tol");
	ADD_PARAMETER(m_abstol      , "abstol");
	ADD_PARAMETER(m_nsolver     , "linear_solver");
	ADD_PARAMETER(m_nschurSolver, "schur_solver");
	ADD_PARAMETER(m_k           , "k");
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
	}

	void* Create(FEModel* fem) override
	{
		FGMRES_Schur_Solver* ls = new FGMRES_Schur_Solver(fem);
		ls->SetPrintLevel(m_print_level);
		ls->SetMaxIterations(m_maxiter);
		ls->SetRelativeResidualTolerance(m_tol);
		ls->ZeroDBlock(m_bzeroDBlock);
		return ls;
	}

private:
	int		m_maxiter;		// max nr of iterations
	int		m_print_level;	// output level
	double	m_tol;
	bool	m_bzeroDBlock;

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FGMRESSchurLinearSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_print_level, "print_level");
	ADD_PARAMETER(m_maxiter    , "maxiter");
	ADD_PARAMETER(m_tol        , "tol");
	ADD_PARAMETER(m_bzeroDBlock, "zero_D_block");
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

	void* Create(FEModel* fem) override
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


} // namespace NumCore

//=============================================================================
// Call this to initialize the NumCore module
void NumCore::InitModule()
{
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

	// register linear solvers
	REGISTER_FECORE_CLASS(SkylineSolver    , "skyline"   );
	REGISTER_FECORE_CLASS(LUSolver         , "LU"        );

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
