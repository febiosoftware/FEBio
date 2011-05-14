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
	FEM& fem = dynamic_cast<FEM&>(*m_pfem);

	// read input data
	if (fem.Input(szfile) == false) return false;

	// initialize and check data
	if (fem.Init() == false) return false;

	// Solve the problem and return error code
	return fem.Solve();
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
bool optimize(FEM& fem, const char* szfile);

bool FEBioOptimize::Run(const char *szfile)
{
	return optimize(dynamic_cast<FEM&>(*m_pfem), szfile);
}

//-----------------------------------------------------------------------------
bool diagnose(FEM& fem, const char* szfile);

bool FEBioDiagnostic::Run(const char *szfile)
{
	return diagnose(dynamic_cast<FEM&>(*m_pfem), szfile);
}
