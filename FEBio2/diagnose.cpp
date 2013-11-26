// diagnose.cpp : implementation of a diagnostic routine
// 
//////////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "FEDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FEBioLib/FEBox.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
//! The diagnose() function performs a diagnostic test on FEBio. 

bool diagnose(FEModel& fem, const char* szfile)
{
	// read the diagnostic file
	// this will also create a specific diagnostic test
	FEDiagnosticImport im;
	FEDiagnostic* pdia = im.LoadFile(fem, szfile);
	if (pdia == 0)
	{
		fprintf(stderr, "Failed reading diagnostic file\n");
		return false;
	}

	// intialize diagnostic
	if (pdia->Init() == false)
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

	// --- run the diagnostic ---

	// the return value will designate the pass/fail result
	bool bret = false;
	try
	{
		bret = pdia->Run();
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
