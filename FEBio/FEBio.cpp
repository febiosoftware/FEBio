///////////////////////////////////////////////////////////////////////////////
//         ________    _________   _________     __     _________            //
//        |        |\ |        |\ |        |\   |  |\  /         \\          //
//        |    ____|| |    ____|| |    __  ||   |__|| |    ___    ||         //
//        |   |\___\| |   |\___\| |   |\_| ||    \_\| |   //  \   ||         //
//        |   ||      |   ||      |   || | ||    __   |  ||    |  ||         //
//        |   ||__    |   ||__    |   ||_| ||   |  |\ |  ||    |  ||         //
//        |       |\  |       |\  |         \\  |  || |  ||    |  ||         //
//        |    ___||  |    ___||  |    ___   || |  || |  ||    |  ||         //
//        |   |\__\|  |   |\__\|  |   |\__|  || |  || |  ||    |  ||         //
//        |   ||      |   ||      |   ||  |  || |  || |  ||    |  ||         //
//        |   ||      |   ||___   |   ||__|  || |  || |   \\__/   ||         //
//        |   ||      |        |\ |          || |  || |           ||         //
//        |___||      |________|| |__________|| |__||  \_________//          //
//                                                                           //
//      F I N I T E   E L E M E N T S   F O R   B I O M E C H A N I C S      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//
//! \mainpage FEBio Documentation
//! \section sec_intro Introduction
//! FEBio is a finite element solver that is specifically designed for three
//! dimensional biomechanical applications. It solves the nonlinear finite
//! element equations using a quasi-Newton method called the BFGS-method. It
//! also offers several biologically relevant constitutive models.
//!
//! This software is developed at the Musculoskeletal Research Laboratories
//! at the University of Utah. All rights reserved.
//! Copyright (c) 2006 - 2010
//!
//! The subversion (svn) revision number of this code can be found in the file
//! FEBio/svnrev.h
//!
//! Main developers:
//!  - Steve Maas
//!  - Jeff Weiss
//!  - Dave Rawlins
//!  - Gerard Ateshian
//!
//! Contributors:
//!  - Alexander Veress
//
///////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "fem.h"
#include "version.h"
#include "FEBioCommand.h"
#include "FECore/FECore.h"
#include "validate.h"
#include "console.h"
#include "log.h"

#define MAXFILE 512

//-----------------------------------------------------------------------------
// FEBio analysis modes
enum FEBIO_MODE {
	FE_UNKNOWN,
	FE_ANALYZE,
	FE_RESTART,
	FE_DIAGNOSE,
	FE_OPTIMIZE
};

//-----------------------------------------------------------------------------
//!  Command line options

//! This structures stores the command line options that were input by the user
struct CMDOPTIONS
{
	FEBIO_MODE	nmode;			//!< FEBio analysis mode
	bool		bdebug;			//!< debug flag

	bool	bsplash;			//!< show splash screen or not
	bool	bsilent;			//!< run FEBio in silent mode (no output to screen)

	char	szfile[MAXFILE];	//!< input file name
	char	szlog [MAXFILE];	//!< log file name
	char	szplt [MAXFILE];	//!< plot file name
	char	szdmp [MAXFILE];	//!< dump file name
	char	szcnf [MAXFILE];	//!< configuration file
};

//-----------------------------------------------------------------------------
// forward declarations
//
bool ParseCmdLine(int argc, char* argv[], CMDOPTIONS& ops);
void Hello(FILE* fp);
int prompt(CMDOPTIONS& ops);

bool optimize(FEM& fem, const char* szfile);
bool diagnose(FEM& fem, const char* szfile);
bool restart (FEM& fem, const char* szfile);
bool solve   (FEM& fem, const char* szfile);

void init_framework(FEM& fem);

int get_app_path (char *pname, size_t pathsize);

//-----------------------------------------------------------------------------
// The starting point of the application
//
int main(int argc, char* argv[])
{
	// parse the command line
	CMDOPTIONS ops;
	if (ParseCmdLine(argc, argv, ops) == false) return 0;

	// load the license file
	LoadLicenseFile();

	// print welcome message
	if (ops.bsplash && (!ops.bsilent)) Hello(stdout);

	// if silent mode only output to file
	if (ops.bsilent) GetLogfile().SetMode(Logfile::FILE_ONLY);

	// if there are no arguments, print the FEBio prompt
	if (argc == 1)
	{
		 if (prompt(ops) == 0) return 0;
	}

	// create the one and only FEM object
	FEM fem;

	// intialize the framework
	init_framework(fem);

	// read the configration file if specified
	if (ops.szcnf[0])
	{
		if (fem.Configure(ops.szcnf) == false) return 1;
	}

	// set the output filenames
	fem.SetLogFilename (ops.szlog);
	fem.SetPlotFilename(ops.szplt);
	fem.SetDumpFilename(ops.szdmp);

	// set options that were passed on the command line
	fem.SetDebugFlag(ops.bdebug);

	// run the FEBio analysis
	bool bret = false;
	switch (ops.nmode)
	{
	case FE_ANALYZE : bret = solve   (fem, ops.szfile); break;
	case FE_OPTIMIZE: bret = optimize(fem, ops.szfile); break;
	case FE_DIAGNOSE: bret = diagnose(fem, ops.szfile); break;
	case FE_RESTART : bret = restart (fem, ops.szfile); break;
	};

	return (bret?0:1);
}

//-----------------------------------------------------------------------------
//!  Parses the command line and returns a CMDOPTIONS structure
//
bool ParseCmdLine(int nargs, char* argv[], CMDOPTIONS& ops)
{
	// set default options
	ops.nmode = FE_UNKNOWN;
	ops.bdebug = false;
	ops.bsplash = true;
	ops.bsilent = false;

	bool blog = false;
	bool bplt = false;
	bool bdmp = false;
	bool brun = true;

	ops.szfile[0] = 0;

	// set the location of the configuration file
	char szpath[1024] = {0};
	get_app_path (szpath, 1023);

	char* ch = strrchr(szpath, '\\');
	if (ch == 0) ch = strrchr(szpath, '/');
	if (ch) ch[1] = 0;

	sprintf(ops.szcnf, "%sfebio.xml", szpath);

	// loop over the arguments
	char* sz;
	for (int i=1; i<nargs; ++i)
	{
		sz = argv[i];

		if (strcmp(sz,"-r") == 0)
		{
			ops.nmode = FE_RESTART;
			strcpy(ops.szfile, argv[++i]);
		}
		else if (strcmp(sz, "-d") == 0)
		{
			ops.nmode = FE_DIAGNOSE;
			strcpy(ops.szfile, argv[++i]);
		}
		else if (strcmp(sz, "-p") == 0)
		{
			bplt = true;
			strcpy(ops.szplt, argv[++i]);
		}
		else if (strcmp(sz, "-a") == 0)
		{
			bdmp = true;
			strcpy(ops.szdmp, argv[++i]);
		}
		else if (strcmp(sz, "-o") == 0)
		{
			blog = true;
			strcpy(ops.szlog, argv[++i]);
		}
		else if (strcmp(sz, "-i") == 0)
		{
			ops.nmode = FE_ANALYZE;
			strcpy(ops.szfile, argv[++i]);
		}
		else if (strcmp(sz, "-s") == 0)
		{
			ops.nmode = FE_OPTIMIZE;
			strcpy(ops.szfile, argv[++i]);
		}
		else if (strcmp(sz, "-g") == 0)
		{
			ops.bdebug = true;
		}
		else if (strcmp(sz, "-nosplash") == 0)
		{
			// don't show the welcome message
			ops.bsplash = false;
		}
		else if (strcmp(sz, "-silent") == 0)
		{
			// no output to screen
			ops.bsilent = true;
		}
		else if (strcmp(sz, "-cnf") == 0)
		{
			strcpy(ops.szcnf, argv[++i]);
		}
		else if (strcmp(sz, "-noconfig") == 0)
		{
			ops.szcnf[0] = 0;
		}
		else if (strcmp(sz, "-info")==0)
		{
			FILE* fp = stdout;
			if ((i<nargs-1) && (argv[i+1][0] != '-'))
			{
				fp = fopen(argv[++i], "wt");
				if (fp == 0) fp = stdout;
			}
			fprintf(fp, "compiled on " __DATE__ "\n");
			fprintf(fp, "FEBio version  = %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
			fprintf(fp, "FECore version = %s\n", FECore::get_version_string());
			fprintf(fp, "SVN revision   = %d\n", SVNREVISION);
			if (fp != stdout) fclose(fp);
		}
		else if (strcmp(sz, "-norun") == 0)
		{
			brun = false;
		}
		else
		{
			fprintf(stderr, "FATAL ERROR: Invalid command line option\n\n");
			return false;
		}
	}

	// derive the other filenames
	char szbase[256]; strcpy(szbase, ops.szfile);
	ch = strrchr(szbase, '.');
	if (ch) *ch = 0;

	char* szext = (ch?ch+1:0);

	strcpy(szpath, ops.szfile);
	ch = strrchr(szpath, '/');
	if (ch == 0) ch = strrchr(szpath, '\\');
	if (ch) *(ch+1) = 0; else szpath[0] = 0;

	if (szext && ((strcmp(szext, "feb")==0) ||
				  (strcmp(szext, "xml")==0) ||
				  (strcmp(szext, "XML")==0)))
	{
		if (!blog) sprintf(ops.szlog, "%s.log", szbase);
		if (!bplt) sprintf(ops.szplt, "%s.plt", szbase);
		if (!bdmp) sprintf(ops.szdmp, "%s.dmp", szbase);
	}
	else
	{
		if (!blog) sprintf(ops.szlog, "%sf3log.txt", szpath);
		if (!bplt) sprintf(ops.szplt, "%sf3plot"   , szpath);
		if (!bdmp) sprintf(ops.szdmp, "%sf3dmp"    , szpath);
	}

	return brun;
}

//-----------------------------------------------------------------------------
//! Initializes the framework
void init_framework(FEM& fem)
{
	FEBioCommand::SetFEM(&fem);
}

//-----------------------------------------------------------------------------
//! Prints the FEBio prompt. If the user did not enter anything on the command
//! line when running FEBio then commands can be entered at the FEBio prompt.
//! This function returns the command arguments as a CMDOPTIONS structure.
int prompt(CMDOPTIONS& ops)
{
	// get a pointer to the console window
	Console* pShell = Console::GetHandle();

	int nargs;
	char* argv[32];

	fprintf(stderr, "Type help for an overview of commands.\n");

	while (1)
	{
		// get a command from the shell
		pShell->GetCommand(nargs, argv);
		if (nargs > 0)
		{
			if (strcmp(argv[0], "quit") == 0) return 0;
			else if (strcmp(argv[0], "help") == 0)
			{
				fprintf(stderr, "\n");
				fprintf(stderr, "help - print this info\n");
				fprintf(stderr, "quit - exits the application\n");
				fprintf(stderr, "run [-i,-s] <file> [OPTIONS] - run an FEBio input file\n");
				fprintf(stderr, "version - print version information\n");
			}
			else if (strcmp(argv[0], "run") == 0)
			{
				ParseCmdLine(nargs, argv, ops);
				return 1;
			}
			else if (strcmp(argv[0], "version") == 0)
			{
				fprintf(stderr, "\nFEBio version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
				fprintf(stderr, "compiled on " __DATE__ "\n");
				fprintf(stderr, "using FECore version %s\n\n", FECore::get_version_string());
			}
			else
			{
				printf("Unknown command: %s\n", argv[0]);
			}
		}
	}
	return 0;
}
