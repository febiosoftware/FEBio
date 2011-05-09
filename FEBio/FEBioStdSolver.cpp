#include "stdafx.h"
#include "FEBioStdSolver.h"
#include "fem.h"

//-----------------------------------------------------------------------------
// This simply calls the FEM::Solve function which will solve the FE problem.
bool FEBioStdSolver::Run(FEModel& mdl, const char* szfile)
{
	FEM& fem = dynamic_cast<FEM&>(mdl);

	// read input data
	if (fem.Input(szfile) == false) return false;

	// initialize and check data
	if (fem.Init() == false) return false;

	// Solve the problem and return error code
	return fem.Solve();
}

//-----------------------------------------------------------------------------
bool FEBioRestart::Run(FEModel& mdl, const char *szfile)
{
	FEM& fem = dynamic_cast<FEM&>(mdl);

	// load restart data
	if (fem.Restart(szfile) == false) return false;

	// continue the analysis
	return fem.Solve();
}

//-----------------------------------------------------------------------------
bool optimize(FEM& fem, const char* szfile);

bool FEBioOptimize::Run(FEModel& fem, const char *szfile)
{
	return optimize(dynamic_cast<FEM&>(fem), szfile);
}

//-----------------------------------------------------------------------------
bool diagnose(FEM& fem, const char* szfile);

bool FEBioDiagnostic::Run(FEModel& fem, const char *szfile)
{
	return diagnose(dynamic_cast<FEM&>(fem), szfile);
}
