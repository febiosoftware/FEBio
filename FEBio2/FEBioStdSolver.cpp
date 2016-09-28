#include "stdafx.h"
#include "FEBioStdSolver.h"
#include "FEBioLib/FEBioModel.h"
#include "FEBioTest/FEDiagnostic.h"
#include "FEBioTest/FETangentDiagnostic.h"
#include "FEBioTest/FERestartDiagnostics.h"
#include "FEBioLib/FEBox.h"
#include "FECore/log.h"
#include "FEBioXML/FERestartImport.h"
#include "FECore/DumpFile.h"

//-----------------------------------------------------------------------------
REGISTER_FECORE_CLASS(FEBioStdSolver     , FETASK_ID, "solve"       );
REGISTER_FECORE_CLASS(FEBioRestart       , FETASK_ID, "restart"     );
REGISTER_FECORE_CLASS(FEBioDiagnostic    , FETASK_ID, "diagnose"    );
REGISTER_FECORE_CLASS(FERestartDiagnostic, FETASK_ID, "restart_test");

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

	// check the extension of the file
	// if the extension is .dmp or not given it is assumed the file
	// is a bindary archive (dump file). Otherwise it is assumed the
	// file is a restart input file.
	const char* ch = strrchr(szfile, '.');
	if ((ch == 0) || (strcmp(ch, ".dmp") == 0) || (strcmp(ch, ".DMP") == 0))
	{
		// the file is binary so just read the dump file and return

		// open the archive
		DumpFile ar(fem);
		if (ar.Open(szfile) == false) { fprintf(stderr, "FATAL ERROR: failed opening restart archive\n"); return false; }

		// read the archive
		try
		{
			if (fem.Serialize(ar) == false) { fprintf(stderr, "FATAL ERROR: failed reading restart data from archive %s\n", szfile); return false; }
		}
		catch (...)
		{
			fprintf(stderr, "FATAL ERROR: failed reading restart data from archive %s\n", szfile); 
			return false;
		}
	}
	else
	{
		// the file is assumed to be a xml-text input file
		FERestartImport file;
		if (file.Load(fem, szfile) == false)
		{
			char szerr[256];
			file.GetErrorMessage(szerr);
			fprintf(stderr, "%s", szerr);
			return false;
		}

		// see if user redefined restart file name
		if (file.m_szdmp[0]) fem.SetDumpFilename(file.m_szdmp);
	}

	// Open the log file for appending
	const char* szlog = fem.GetLogfileName();
	if (felog.append(szlog) == false)
	{
		printf("WARNING: Could not reopen log file. A new log file is created\n");
		felog.open(szlog);
		return false;
	}

	// inform the user from where the problem is restarted
	felog.printbox(" - R E S T A R T -", "Restarting from time %lg.\n", fem.m_ftime);

	return true;
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
		felog.SetMode(Logfile::LOG_FILE_AND_SCREEN);
		felog.printf("Exception thrown. Aborting diagnostic.\n");
		bret = false;
	}

	if (bret) felog.printf("Diagnostic passed\n");
	else felog.printf("Diagnostic failed\n");

	return bret;
}
