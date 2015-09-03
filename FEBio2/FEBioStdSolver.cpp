#include "stdafx.h"
#include "FEBioStdSolver.h"
#include "FEBioLib/FEBioModel.h"
#include "FEBioTest/FEDiagnostic.h"
#include "FEBioTest/FETangentDiagnostic.h"
#include "FEBioLib/FEBox.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
REGISTER_FECORE_CLASS(FEBioStdSolver , FETASK_ID, "solve"   );
REGISTER_FECORE_CLASS(FEBioRestart   , FETASK_ID, "restart" );
REGISTER_FECORE_CLASS(FEBioDiagnostic, FETASK_ID, "diagnose");

//-----------------------------------------------------------------------------
FEBioStdSolver::FEBioStdSolver(FEModel* pfem) : FECoreTask(pfem) {}

//-----------------------------------------------------------------------------
// This simply calls the FEModel::Init
bool FEBioStdSolver::Init(const char* szfile)
{
	return (m_pfem ? m_pfem->Init() : false);
}

//-----------------------------------------------------------------------------
// This simply calls the FEM::Solve function which will solve the FE problem.
bool FEBioStdSolver::Run()
{
	// Solve the problem and return error code
	return (m_pfem ? m_pfem->Solve() : false);
}

//-----------------------------------------------------------------------------
bool FEBioRestart::Init(const char *szfile)
{
	FEBioModel& fem = static_cast<FEBioModel&>(*GetFEModel());

	// load restart data
	return fem.Restart(szfile);
}

//-----------------------------------------------------------------------------
bool FEBioRestart::Run()
{
	// continue the analysis
	return (m_pfem ? m_pfem->Solve() : false);
}

//-----------------------------------------------------------------------------
bool FEBioDiagnostic::Init(const char *szfile)
{
	FEModel& fem = *GetFEModel();

	// read the diagnostic file
	// this will also create a specific diagnostic test
	FEDiagnosticImport im;
	m_pdia = im.LoadFile(fem, szfile);
	if (m_pdia == 0)
	{
		fprintf(stderr, "Failed reading diagnostic file\n");
		return false;
	}

	// intialize diagnostic
	if (m_pdia->Init() == false)
	{
		fprintf(stderr, "Diagnostic initialization failed\n\n");
		return false;
	}

	// --- initialize FE Model data ---
	if (fem.Init() == false)
	{
		fprintf(stderr, "FE-model data initialized has failed\n\n");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEBioDiagnostic::Run()
{
	if (m_pdia == 0) return false;
	
	// --- run the diagnostic ---

	// the return value will designate the pass/fail result
	bool bret = false;
	try
	{
		bret = m_pdia->Run();
	}
	catch (...)
	{
		felog.SetMode(Logfile::FILE_AND_SCREEN);
		felog.printf("Exception thrown. Aborting diagnostic.\n");
		bret = false;
	}

	if (bret) felog.printf("Diagnostic passed\n");
	else felog.printf("Diagnostic failed\n");

	return bret;
}
