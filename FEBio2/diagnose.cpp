// diagnose.cpp : implementation of a diagnostic routine
// 
//////////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "FEDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FEBox.h"
#include "FEBioLib/log.h"

//-----------------------------------------------------------------------------
//! The diagnose() function performs a diagnostic test on FEBio. 

bool diagnose(FEM& fem, const char* szfile)
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

	// --- initialize FEM data ---
	if (fem.Init() == false)
	{
		fprintf(stderr, "FEM initialized has failed\n\n");
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
		clog.SetMode(Logfile::FILE_AND_SCREEN);
		clog.printf("Exception thrown. Aborting diagnostic.\n");
		bret = false;
	}

	if (bret) clog.printf("Diagnostic passed\n");
	else clog.printf("Diagnostic failed\n");

	return bret;
}
