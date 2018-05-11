#include "stdafx.h"
#include "NumCore.h"
#include "SkylineSolver.h"
#include "LUSolver.h"
#include "ConjGradIterSolver.h"
#include "PSLDLTSolver.h"
#include "SuperLUSolver.h"
#include "SuperLU_MT_Solver.h"
#include "PardisoSolver.h"
#include "WSMPSolver.h"
#include "RCICGSolver.h"
#include "FGMRESSolver.h"
#include "FGMRES_ILU0_Solver.h"
#include "FGMRES_ILUT_Solver.h"
#include "BIPNSolver.h"
#include "HypreGMRESsolver.h"
#include "FECore/FE_enum.h"
#include "FECore/FECoreFactory.h"
#include "FECore/FECoreKernel.h"

namespace NumCore {

template <class T, int nid> class LinearSolverFactory_T : public FELinearSolverFactory
{
public:
	LinearSolverFactory_T() : FELinearSolverFactory(nid)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterLinearSolver(this);
	}
	LinearSolver* Create() { return new T(); }
};

template <> class LinearSolverFactory_T<RCICGSolver, RCICG_SOLVER> : public FELinearSolverFactory
{
public:
	LinearSolverFactory_T() : FELinearSolverFactory(RCICG_SOLVER)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterLinearSolver(this);

		m_maxiter = 0;
		m_tol = 1e-5;
		m_precond = 1;
		m_print_level = 0;
	}

	LinearSolver* Create() override
	{
		RCICGSolver* ls = new RCICGSolver();
		ls->SetMaxIterations(m_maxiter);
		ls->SetTolerance(m_tol);
		ls->SetPreconditioner(m_precond);
		ls->SetPrintLevel(m_print_level);
		return ls;
	}

private:
	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_precond;		// pre-conditioner
	int		m_print_level;	// output level

	DECLARE_PARAMETER_LIST();
};

typedef LinearSolverFactory_T<RCICGSolver, RCICG_SOLVER> RCICG_SolverFactory;

BEGIN_PARAMETER_LIST(RCICG_SolverFactory, FELinearSolverFactory)
	ADD_PARAMETER(m_maxiter, FE_PARAM_INT, "maxiter");
	ADD_PARAMETER(m_tol, FE_PARAM_DOUBLE, "tol");
	ADD_PARAMETER(m_precond, FE_PARAM_INT, "precondition");
	ADD_PARAMETER(m_print_level, FE_PARAM_INT, "print_level");
END_PARAMETER_LIST();

template <> class LinearSolverFactory_T<FGMRES_ILUT_Solver, FGMRES_ILUT_SOLVER> : public FELinearSolverFactory
{
public:
	LinearSolverFactory_T() : FELinearSolverFactory(FGMRES_ILUT_SOLVER)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterLinearSolver(this);

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
	LinearSolver* Create() override
	{ 
		FGMRES_ILUT_Solver* ls = new FGMRES_ILUT_Solver();
		ls->SetMaxFill(m_maxfill);
		ls->SetFillTolerance(m_fillTol);
		ls->SetMaxIterations(m_maxiter);
		ls->SetNonRestartedIterations(m_nrestart);
		ls->SetPrintLevel(m_print_level);
		ls->DoResidualStoppingTest(m_doResidualTest);
		ls->SetResidualTolerance(m_tol);

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

	DECLARE_PARAMETER_LIST();
};

typedef LinearSolverFactory_T<FGMRES_ILUT_Solver, FGMRES_ILUT_SOLVER> FGMRES_ILUT_SolverFactory;

BEGIN_PARAMETER_LIST(FGMRES_ILUT_SolverFactory, FELinearSolverFactory)
	ADD_PARAMETER(m_maxfill, FE_PARAM_INT, "maxfill");
	ADD_PARAMETER(m_fillTol, FE_PARAM_DOUBLE, "filltol");
	ADD_PARAMETER(m_maxiter       , FE_PARAM_INT   , "maxiter");
	ADD_PARAMETER(m_nrestart      , FE_PARAM_INT   , "maxrestart");
	ADD_PARAMETER(m_print_level   , FE_PARAM_INT   , "print_level");
	ADD_PARAMETER(m_doResidualTest, FE_PARAM_BOOL  , "check_residual");
	ADD_PARAMETER(m_tol           , FE_PARAM_DOUBLE, "tol");
	ADD_PARAMETER(m_checkZeroDiagonal, FE_PARAM_BOOL  , "replace_zero_diagonal");
	ADD_PARAMETER(m_zeroThreshold    , FE_PARAM_DOUBLE, "zero_threshold");
	ADD_PARAMETER(m_zeroReplace      , FE_PARAM_DOUBLE, "zero_replace");
END_PARAMETER_LIST();

template <> class LinearSolverFactory_T<FGMRES_ILU0_Solver, FGMRES_ILU0_SOLVER> : public FELinearSolverFactory
{
public:
	LinearSolverFactory_T() : FELinearSolverFactory(FGMRES_ILU0_SOLVER)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterLinearSolver(this);

		m_maxiter = 0; // use default min(N, 150)
		m_nrestart = 0;
		m_print_level = 0;
		m_doResidualTest = true;
		m_tol = 0;

		m_checkZeroDiagonal = true;
		m_zeroThreshold = 1e-16;
		m_zeroReplace = 1e-10;
	}
	LinearSolver* Create() override
	{
		FGMRES_ILU0_Solver* ls = new FGMRES_ILU0_Solver();
		ls->SetMaxIterations(m_maxiter);
		ls->SetNonRestartedIterations(m_nrestart);
		ls->SetPrintLevel(m_print_level);
		ls->DoResidualStoppingTest(m_doResidualTest);
		ls->SetResidualTolerance(m_tol);

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

	// pre-conditioner parameters
	bool	m_checkZeroDiagonal;	// check for zero diagonals
	double	m_zeroThreshold;		// threshold for zero diagonal check
	double	m_zeroReplace;			// replacement value for zero diagonal

	DECLARE_PARAMETER_LIST();
};

typedef LinearSolverFactory_T<FGMRES_ILU0_Solver, FGMRES_ILU0_SOLVER> FGMRES_ILU0_SolverFactory;

BEGIN_PARAMETER_LIST(FGMRES_ILU0_SolverFactory, FELinearSolverFactory)
	ADD_PARAMETER(m_maxiter       , FE_PARAM_INT   , "maxiter");
	ADD_PARAMETER(m_nrestart      , FE_PARAM_INT   , "maxrestart");
	ADD_PARAMETER(m_print_level   , FE_PARAM_INT   , "print_level");
	ADD_PARAMETER(m_doResidualTest, FE_PARAM_BOOL  , "check_residual");
	ADD_PARAMETER(m_tol           , FE_PARAM_DOUBLE, "tol");
	ADD_PARAMETER(m_checkZeroDiagonal, FE_PARAM_BOOL  , "replace_zero_diagonal");
	ADD_PARAMETER(m_zeroThreshold    , FE_PARAM_DOUBLE, "zero_threshold");
	ADD_PARAMETER(m_zeroReplace      , FE_PARAM_DOUBLE, "zero_replace");
END_PARAMETER_LIST();


template <> class LinearSolverFactory_T<FGMRESSolver, FGMRES_SOLVER> : public FELinearSolverFactory
{
public:
	LinearSolverFactory_T() : FELinearSolverFactory(FGMRES_SOLVER)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterLinearSolver(this);

		m_maxiter = 0; // use default min(N, 150)
		m_nrestart = 0; // use maxiter
		m_print_level = 0;
		m_doResidualTest = true;
		m_tol = 0;
	}
	LinearSolver* Create() override
	{
		FGMRESSolver* ls = new FGMRESSolver();
		ls->SetMaxIterations(m_maxiter);
		ls->SetNonRestartedIterations(m_nrestart);
		ls->SetPrintLevel(m_print_level);
		ls->DoResidualStoppingTest(m_doResidualTest);
		ls->SetResidualTolerance(m_tol);
		return ls;
	}

private:
	int		m_maxiter;			// max number of iterations
	int		m_nrestart;			// nr of non-restarted iterations
	int		m_print_level;		// print level
	bool	m_doResidualTest;	// residual stopping tets flag
	double	m_tol;				// residual convergence tolerance

	DECLARE_PARAMETER_LIST();
};

typedef LinearSolverFactory_T<FGMRESSolver, FGMRES_SOLVER> FGMRESSolverFactory;

BEGIN_PARAMETER_LIST(FGMRESSolverFactory, FELinearSolverFactory)
	ADD_PARAMETER(m_maxiter       , FE_PARAM_INT   , "maxiter");
	ADD_PARAMETER(m_print_level   , FE_PARAM_INT   , "print_level");
	ADD_PARAMETER(m_doResidualTest, FE_PARAM_BOOL  , "check_residual");
	ADD_PARAMETER(m_nrestart      , FE_PARAM_INT   , "maxrestart");
	ADD_PARAMETER(m_tol           , FE_PARAM_DOUBLE, "tol");
END_PARAMETER_LIST();

#define REGISTER_LINEAR_SOLVER(theSolver, theID) static LinearSolverFactory_T<theSolver, theID> _##theSolver;

//=======================================================================================

template <> class LinearSolverFactory_T<BIPNSolver, BIPN_SOLVER> : public FELinearSolverFactory
{
public:
	LinearSolverFactory_T() : FELinearSolverFactory(BIPN_SOLVER)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterLinearSolver(this);

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

	LinearSolver* Create() override
	{
		BIPNSolver* ls = new BIPNSolver();
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

	DECLARE_PARAMETER_LIST();
};

typedef LinearSolverFactory_T<BIPNSolver, BIPN_SOLVER> BIPN_SolverFactory;

BEGIN_PARAMETER_LIST(BIPN_SolverFactory, FELinearSolverFactory)
	ADD_PARAMETER(m_maxiter    , FE_PARAM_INT   , "maxiter"    );
	ADD_PARAMETER(m_tol        , FE_PARAM_DOUBLE, "tol"        );
	ADD_PARAMETER(m_print_level, FE_PARAM_INT   , "print_level");
	ADD_PARAMETER(m_use_cg     , FE_PARAM_BOOL  , "use_cg");
	ADD_PARAMETER(m_cg_max     , FE_PARAM_INT   , "cg_maxiter" );
	ADD_PARAMETER(m_cg_tol     , FE_PARAM_DOUBLE, "cg_tol"     );
	ADD_PARAMETER(m_cg_res     , FE_PARAM_BOOL  , "cg_check_residual");
	ADD_PARAMETER(m_gmres_max  , FE_PARAM_INT   , "gmres_maxiter");
	ADD_PARAMETER(m_gmres_tol  , FE_PARAM_DOUBLE, "gmres_tol" );
	ADD_PARAMETER(m_gmres_res  , FE_PARAM_BOOL  , "gmres_check_residual");
	ADD_PARAMETER(m_gmres_ilu0 , FE_PARAM_BOOL  , "gmres_precondition");
END_PARAMETER_LIST();

//=======================================================================================

template <> class LinearSolverFactory_T<HypreGMRESsolver, HYPRE_GMRES> : public FELinearSolverFactory
{
public:
	LinearSolverFactory_T() : FELinearSolverFactory(HYPRE_GMRES)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterLinearSolver(this);

		m_maxiter = 1000;
		m_tol = 1e-7;
		m_print_level = 0;
	}

	LinearSolver* Create() override
	{
		HypreGMRESsolver* ls = new HypreGMRESsolver();
		ls->SetPrintLevel(m_print_level);
		ls->SetMaxIterations(m_maxiter);
		ls->SetConvergencTolerance(m_tol);
		return ls;
	}

private:
	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_print_level;	// output level

	DECLARE_PARAMETER_LIST();
};

typedef LinearSolverFactory_T<HypreGMRESsolver, HYPRE_GMRES> HypreGMRES_SolverFactory;

BEGIN_PARAMETER_LIST(HypreGMRES_SolverFactory, FELinearSolverFactory)
	ADD_PARAMETER(m_print_level, FE_PARAM_INT, "print_level");
	ADD_PARAMETER(m_maxiter    , FE_PARAM_INT   , "maxiter"    );
	ADD_PARAMETER(m_tol        , FE_PARAM_DOUBLE, "tol"        );
END_PARAMETER_LIST();


} // namespace NumCore

//=============================================================================
// Call this to initialize the NumCore module
void NumCore::InitModule()
{
REGISTER_LINEAR_SOLVER(SkylineSolver     , SKYLINE_SOLVER     );
REGISTER_LINEAR_SOLVER(PSLDLTSolver      , PSLDLT_SOLVER      );
REGISTER_LINEAR_SOLVER(SuperLUSolver     , SUPERLU_SOLVER     );
REGISTER_LINEAR_SOLVER(SuperLU_MT_Solver , SUPERLU_MT_SOLVER  );
REGISTER_LINEAR_SOLVER(PardisoSolver     , PARDISO_SOLVER     );
REGISTER_LINEAR_SOLVER(LUSolver          , LU_SOLVER          );
REGISTER_LINEAR_SOLVER(WSMPSolver        , WSMP_SOLVER        );
REGISTER_LINEAR_SOLVER(ConjGradIterSolver, CG_ITERATIVE_SOLVER);
REGISTER_LINEAR_SOLVER(RCICGSolver       , RCICG_SOLVER       );
REGISTER_LINEAR_SOLVER(FGMRESSolver      , FGMRES_SOLVER      );
REGISTER_LINEAR_SOLVER(FGMRES_ILUT_Solver, FGMRES_ILUT_SOLVER );
REGISTER_LINEAR_SOLVER(FGMRES_ILU0_Solver, FGMRES_ILU0_SOLVER );
REGISTER_LINEAR_SOLVER(BIPNSolver        , BIPN_SOLVER        );
REGISTER_LINEAR_SOLVER(HypreGMRESsolver  , HYPRE_GMRES        );
}
