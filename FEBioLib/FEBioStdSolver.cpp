#include "stdafx.h"
#include "FEBioStdSolver.h"
#include <FEBioLib/FEBioModel.h>
#include <FECore/log.h>
#include <FEBioXML/FERestartImport.h>
#include <FECore/DumpFile.h>

//-----------------------------------------------------------------------------
FEBioStdSolver::FEBioStdSolver(FEModel* pfem) : FECoreTask(pfem) {}

//-----------------------------------------------------------------------------
// This simply calls the FEModel::Init
bool FEBioStdSolver::Init(const char* szfile)
{
	return (GetFEModel() ? GetFEModel()->Init() : false);
}

//-----------------------------------------------------------------------------
// This simply calls the FEM::Solve function which will solve the FE problem.
bool FEBioStdSolver::Run()
{
	// Solve the problem and return error code
	return (GetFEModel() ? GetFEModel()->Solve() : false);
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
			fem.Serialize(ar);
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
	felog.printbox(" - R E S T A R T -", "Restarting from time %lg.\n", fem.GetCurrentTime());

	return true;
}

//-----------------------------------------------------------------------------
bool FEBioRestart::Run()
{
	// continue the analysis
	return (GetFEModel() ? GetFEModel()->Solve() : false);
}
