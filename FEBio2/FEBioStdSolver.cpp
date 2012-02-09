#include "stdafx.h"
#include "FEBioStdSolver.h"
#include "fem.h"

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEBioStdSolver , FEBioTask, "solve"   );
REGISTER_FEBIO_CLASS(FEBioRestart   , FEBioTask, "restart" );
REGISTER_FEBIO_CLASS(FEBioOptimize  , FEBioTask, "optimize");
REGISTER_FEBIO_CLASS(FEBioDiagnostic, FEBioTask, "diagnose");

//-----------------------------------------------------------------------------
// This simply calls the FEM::Solve function which will solve the FE problem.
bool FEBioStdSolver::Run(const char* szfile)
{
	// read input data
	if (m_pfem->Input(szfile) == false) return false;

	// initialize and check data
	if (m_pfem->Init() == false) return false;

	// Solve the problem and return error code
	return m_pfem->Solve();
}

//-----------------------------------------------------------------------------
bool FEBioRestart::Run(const char *szfile)
{
	FEM& fem = dynamic_cast<FEM&>(*m_pfem);

	// load restart data
	if (fem.Restart(szfile) == false) return false;

	// continue the analysis
	return fem.Solve();
}

//-----------------------------------------------------------------------------
bool optimize(FEModel& fem, const char* szfile);

bool FEBioOptimize::Run(const char *szfile)
{
	return optimize(*m_pfem, szfile);
}

//-----------------------------------------------------------------------------
bool diagnose(FEModel& fem, const char* szfile);

bool FEBioDiagnostic::Run(const char *szfile)
{
	return diagnose(*m_pfem, szfile);
}
