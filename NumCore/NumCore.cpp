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
}
