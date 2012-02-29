// restart module
#include "stdafx.h"
#include "fem.h"
#include "FEBioProgress.h"
#include "FEBioXML/FERestartImport.h"
#include "FEBioLib/log.h"

//-----------------------------------------------------------------------------
bool restart(FEM& fem, const char* szfile)
{
	// load restart data
	if (fem.Restart(szfile) == false) return false;

	// create progress tracker
	FEBioProgress prg(fem);

	// continue the analysis
	return fem.Solve(prg);
}

//-----------------------------------------------------------------------------
//!  This routine reads a binary archive that stores a restart point and prepares
//!  the FEM data to be restarted from this point
//!	\param[in] szfile name of the file

bool FEM::Restart(const char* szfile)
{
	// check the extension of the file
	// if the extension is .dmp or not given it is assumed the file
	// is a bindary archive (dump file). Otherwise it is assumed the
	// file is a restart input file.
	const char* ch = strrchr(szfile, '.');
	if ((ch == 0) || (strcmp(ch, ".dmp") == 0) || (strcmp(ch, ".DMP") == 0))
	{
		// the file is binary so just read the dump file and return

		// open the archive
		DumpFile ar(this);
		if (ar.Open(szfile) == false) { fprintf(stderr, "FATAL ERROR: failed opening restart archive\n"); return false; }

		// read the archive
		try
		{
			if (Serialize(ar) == false) { fprintf(stderr, "FATAL ERROR: failed reading restart data from archive %s\n", szfile); return false; }
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
		if (file.Load(*this, szfile) == false)
		{
			char szerr[256];
			file.GetErrorMessage(szerr);
			fprintf(stderr, "%s", szerr);
			return false;
		}

		// see if user redefined restart file name
		if (file.m_szdmp[0]) SetDumpFilename(file.m_szdmp);
	}

	// Open the log file for appending
	if (clog.append(m_szlog) == false)
	{
		printf("WARNING: Could not reopen log file. A new log file is created\n");
		clog.open(m_szlog);
		return false;
	}

	// inform the user from where the problem is restarted
	clog.printbox(" - R E S T A R T -", "Restarting from time %lg.\n", m_ftime);

	return true;
}
