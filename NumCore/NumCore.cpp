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

LinearSolver* NumCore::CreateLinearSolver(int ntype)
{
	LinearSolver* pls = 0;
	switch (ntype)
	{
	case SKYLINE_SOLVER      : pls = new SkylineSolver(); break;
	case PSLDLT_SOLVER       : pls = new PSLDLTSolver (); break;
	case SUPERLU_SOLVER      : pls = new SuperLUSolver(); break;
	case SUPERLU_MT_SOLVER   : pls = new SuperLU_MT_Solver(); break;
	case PARDISO_SOLVER      : pls = new PardisoSolver(); break;
	case LU_SOLVER           : pls = new LUSolver(); break;
	case WSMP_SOLVER         : pls = new WSMPSolver(); break;
	case CG_ITERATIVE_SOLVER : pls = new ConjGradIterSolver(); break;
	case RCICG_SOLVER        : pls = new RCICGSolver(); break;
	}

	return pls;
}
