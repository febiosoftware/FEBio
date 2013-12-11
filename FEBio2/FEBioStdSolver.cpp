#include "stdafx.h"
#include "FEBioStdSolver.h"
#include "FEBioLib/FEBioModel.h"

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEBioStdSolver , FETASK_ID, "solve"   );
REGISTER_FEBIO_CLASS(FEBioRestart   , FETASK_ID, "restart" );
REGISTER_FEBIO_CLASS(FEBioDiagnostic, FETASK_ID, "diagnose");

//-----------------------------------------------------------------------------
// This simply calls the FEM::Solve function which will solve the FE problem.
bool FEBioStdSolver::Run(const char* szfile)
{
	FEBioModel& fem = dynamic_cast<FEBioModel&>(*m_pfem);

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
	FEBioModel& fem = dynamic_cast<FEBioModel&>(*m_pfem);

	// load restart data
	if (fem.Restart(szfile) == false) return false;

	// continue the analysis
	return fem.Solve();
}


//-----------------------------------------------------------------------------
bool diagnose(FEModel& fem, const char* szfile);

bool FEBioDiagnostic::Run(const char *szfile)
{
	return diagnose(*m_pfem, szfile);
}
