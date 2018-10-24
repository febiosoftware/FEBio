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
#include "StokesSolver.h"
#include "CG_Stokes_Solver.h"
#include "SchurSolver.h"
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
class CGStokesSolverFactory : public LinearSolverFactory
{
public:
	CGStokesSolverFactory() : LinearSolverFactory("cg_stokes")
	{
		m_maxiter = 0;
		m_tol = 1e-5;
		m_print_level = 0;
	}

	void* Create(FEModel* fem) override
	{
		CG_Stokes_Solver* ls = new CG_Stokes_Solver(fem);
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

BEGIN_FECORE_CLASS(CGStokesSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_maxiter,  "maxiter");
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

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FGMRES_ILU0_Factory, LinearSolverFactory)
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
		ls->SetResidualTolerance(m_tol);
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
	HYPRE_FGMRES_SolverFactory() : LinearSolverFactory("hypre_fgmres")
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

class StokesLinearSolverFactory : public LinearSolverFactory
{
public:
	StokesLinearSolverFactory() : LinearSolverFactory("stokes")
	{
		m_maxiter = 0;
		m_tol = 1e-7;
		m_print_level = 0;
	}

	void* Create(FEModel* fem) override
	{
		StokesSolver* ls = new StokesSolver(fem);
		ls->SetPrintLevel(m_print_level);
		ls->SetMaxIterations(m_maxiter);
		ls->SetConvergenceTolerance(m_tol);
		return ls;
	}

private:
	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_print_level;	// output level

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(StokesLinearSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_print_level, "print_level");
	ADD_PARAMETER(m_maxiter    , "maxiter");
	ADD_PARAMETER(m_tol        , "tol");
END_FECORE_CLASS();

//=============================================================================

class SchurLinearSolverFactory : public LinearSolverFactory
{
public:
	SchurLinearSolverFactory() : LinearSolverFactory("schur")
	{
		m_maxiter = 0;
		m_tol = 1e-7;
		m_print_level = 0;
	}

	void* Create(FEModel* fem) override
	{
		SchurSolver* ls = new SchurSolver(fem);
		ls->SetPrintLevel(m_print_level);
		ls->SetMaxIterations(m_maxiter);
		ls->SetConvergenceTolerance(m_tol);
		return ls;
	}

private:
	int		m_maxiter;		// max nr of iterations
	double	m_tol;			// residual relative tolerance
	int		m_print_level;	// output level

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(SchurLinearSolverFactory, LinearSolverFactory)
	ADD_PARAMETER(m_print_level, "print_level");
	ADD_PARAMETER(m_maxiter    , "maxiter");
	ADD_PARAMETER(m_tol        , "tol");
END_FECORE_CLASS();


} // namespace NumCore

//=============================================================================
// Call this to initialize the NumCore module
void NumCore::InitModule()
{
	REGISTER_FECORE_FACTORY(SchurLinearSolverFactory  );
	REGISTER_FECORE_FACTORY(StokesLinearSolverFactory );
	REGISTER_FECORE_FACTORY(HYPRE_FGMRES_SolverFactory);
	REGISTER_FECORE_FACTORY(BIPNSolverFactory         );
	REGISTER_FECORE_FACTORY(FGMRESSolverFactory       );
	REGISTER_FECORE_FACTORY(FGMRES_ILU0_Factory       );
	REGISTER_FECORE_FACTORY(FGMRES_ILUT_Factory       );
	REGISTER_FECORE_FACTORY(CGStokesSolverFactory     );
	REGISTER_FECORE_FACTORY(RCICGSolverFactory        );

	REGISTER_FECORE_CLASS(WSMPSolver       , FELINEARSOLVER_ID, "wsmp"      );
	REGISTER_FECORE_CLASS(SuperLUSolver    , FELINEARSOLVER_ID, "superlu"   );
	REGISTER_FECORE_CLASS(SuperLU_MT_Solver, FELINEARSOLVER_ID, "superlu_mt");
	REGISTER_FECORE_CLASS(SkylineSolver    , FELINEARSOLVER_ID, "skyline"   );
	REGISTER_FECORE_CLASS(PSLDLTSolver     , FELINEARSOLVER_ID, "psldlt"    );
	REGISTER_FECORE_CLASS(PardisoSolver    , FELINEARSOLVER_ID, "pardiso"   );
	REGISTER_FECORE_CLASS(LUSolver         , FELINEARSOLVER_ID, "LU"        );

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

