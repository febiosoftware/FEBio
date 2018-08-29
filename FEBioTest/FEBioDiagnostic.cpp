#include "FEBioDiagnostic.h"
#include <FECore/log.h>

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
