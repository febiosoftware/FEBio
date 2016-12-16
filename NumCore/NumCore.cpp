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

template <> class LinearSolverFactory_T<FGMRES_ILUT_Solver, FGMRES_ILUT_SOLVER> : public FELinearSolverFactory
{
public:
	LinearSolverFactory_T() : FELinearSolverFactory(FGMRES_ILUT_SOLVER)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterLinearSolver(this);

		m_fillTol = 1e-6;
		m_maxfill = 1;
		m_maxiter = 0; // use default min(N, 150)
	}
	LinearSolver* Create() 
	{ 
		FGMRES_ILUT_Solver* ls = new FGMRES_ILUT_Solver();
		ls->m_maxfill = m_maxfill;
		ls->m_fillTol = m_fillTol;
		ls->m_maxiter = m_maxiter;
		return ls;
	}

private:
	int		m_maxfill;		// max fill in values (I think this is in terms of bandwidth, not actual values)
	double	m_fillTol;		// tolerance for fill in criterion
	int		m_maxiter;		// max number of iterations

	DECLARE_PARAMETER_LIST();
};

typedef LinearSolverFactory_T<FGMRES_ILUT_Solver, FGMRES_ILUT_SOLVER> FGMRES_ILUT_SolverFactory;

BEGIN_PARAMETER_LIST(FGMRES_ILUT_SolverFactory, FELinearSolverFactory)
	ADD_PARAMETER(m_maxfill, FE_PARAM_INT, "maxfil");
	ADD_PARAMETER(m_fillTol, FE_PARAM_DOUBLE, "tol");
	ADD_PARAMETER(m_maxiter, FE_PARAM_INT, "maxiter");
END_PARAMETER_LIST();

template <> class LinearSolverFactory_T<FGMRES_ILU0_Solver, FGMRES_ILU0_SOLVER> : public FELinearSolverFactory
{
public:
	LinearSolverFactory_T() : FELinearSolverFactory(FGMRES_ILU0_SOLVER)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterLinearSolver(this);

		m_maxiter = 0; // use default min(N, 150)
	}
	LinearSolver* Create()
	{
		FGMRES_ILU0_Solver* ls = new FGMRES_ILU0_Solver();
		ls->m_maxiter = m_maxiter;
		return ls;
	}

private:
	int		m_maxiter;		// max number of iterations

	DECLARE_PARAMETER_LIST();
};

typedef LinearSolverFactory_T<FGMRES_ILU0_Solver, FGMRES_ILU0_SOLVER> FGMRES_ILU0_SolverFactory;

BEGIN_PARAMETER_LIST(FGMRES_ILU0_SolverFactory, FELinearSolverFactory)
	ADD_PARAMETER(m_maxiter, FE_PARAM_INT, "maxiter");
END_PARAMETER_LIST();


template <> class LinearSolverFactory_T<FGMRESSolver, FGMRES_SOLVER> : public FELinearSolverFactory
{
public:
	LinearSolverFactory_T() : FELinearSolverFactory(FGMRES_SOLVER)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		fecore.RegisterLinearSolver(this);

		m_maxiter = 0; // use default min(N, 150)
	}
	LinearSolver* Create()
	{
		FGMRESSolver* ls = new FGMRESSolver();
		ls->m_maxiter = m_maxiter;
		return ls;
	}

private:
	int		m_maxiter;		// max number of iterations

	DECLARE_PARAMETER_LIST();
};

typedef LinearSolverFactory_T<FGMRESSolver, FGMRES_SOLVER> FGMRESSolverFactory;

BEGIN_PARAMETER_LIST(FGMRESSolverFactory, FELinearSolverFactory)
	ADD_PARAMETER(m_maxiter, FE_PARAM_INT, "maxiter");
END_PARAMETER_LIST();

#define REGISTER_LINEAR_SOLVER(theSolver, theID) static LinearSolverFactory_T<theSolver, theID> _##theSolver;
}

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
REGISTER_LINEAR_SOLVER(FGMRES_ILUT_Solver , FGMRES_ILUT_SOLVER  );
REGISTER_LINEAR_SOLVER(FGMRES_ILU0_Solver , FGMRES_ILU0_SOLVER  );
}
