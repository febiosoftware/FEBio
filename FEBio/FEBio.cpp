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
#include "FEBioCommand.h"
#include "FECore/FECore.h"
#include "validate.h"
#include "console.h"

#define MAXFILE 256

//-----------------------------------------------------------------------------
//!  Command line options

//! This structures stores the command line options that were input by the user

struct CMDOPTIONS
{
	char	szfile[MAXFILE];	//!< input file name

	bool	blog;
	char	szlog[MAXFILE];		//!< log file name

	bool	bplt;
	char	szplt[MAXFILE];		//!< plot file name

	bool	bdmp;				//!< dump flag
	char	szdmp[MAXFILE];		//!< dump file name

	char	szcnf[MAXFILE];		//!< configuration file

	bool	bdiag;				//!< run diagnostic file
	bool	brun;				//!< run the problem or just check data?
	bool	brstrt;				//!< restart flag
	bool	boptim;				//!< optimization flag
	bool	bdebug;				//!< debug flag

	bool	bsplash;			//!< show splash screen or not
};

///////////////////////////////////////////////////////////////////////////////
// forward declarations
//

bool ParseCmdLine(int argc, char* argv[], CMDOPTIONS& ops);
void Hello(FILE* fp);
void print_banner();

bool optimize(FEM& fem, const char* szfile);

bool diagnose(FEM& fem, const char* szfile);

void init_framework(FEM& fem);

int get_app_path (char *pname, size_t pathsize);

int prompt(CMDOPTIONS& ops);

//-----------------------------------------------------------------------------
// The starting point of the application

int main(int argc, char* argv[])
{
	// parse the command line
	CMDOPTIONS ops;
	if (ParseCmdLine(argc, argv, ops) == false) return 0;

	// load the license file
	LoadLicenseFile();

	// print welcome message
	if (ops.bsplash)
	{
//#ifdef WIN32
//		print_banner();
//#else
		Hello(stdout);
//#endif
	}

	// if there are no arguments, print the FEBio prompt
	if (argc == 1)
	{
		 int nret = prompt(ops);
		 if (nret == 0) return 0;
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

	// set the filenames
	fem.SetLogFilename (ops.szlog);
	fem.SetPlotFilename(ops.szplt);
	fem.SetDumpFilename(ops.szdmp);

	// set options that were passed on the command line
	fem.SetDebugFlag(ops.bdebug);

	// set the default plot and print levels
/*
	if (ops.bdebug)
	{
		fem.m_pStep->SetPrintLevel(FE_PRINT_MINOR_ITRS_EXP);
		fem.m_pStep->SetPlotLevel (FE_PLOT_MINOR_ITRS);
	}
	else
	{
		fem.m_pStep->SetPrintLevel(FE_PRINT_MINOR_ITRS);
		fem.m_pStep->SetPlotLevel (FE_PLOT_MAJOR_ITRS);
	}
*/
	// check for parameter optimization
	if (ops.boptim) optimize(fem, ops.szfile);
	else if (ops.bdiag) diagnose(fem, ops.szfile);
	else
	{
		// input data
		if (ops.brstrt)
		{
			// do a restart
			if (fem.Restart(ops.szfile) == false) return 1;
		}
		else
		{
			// read input data
			if (fem.Input(ops.szfile) == false) return 1;

			// initialize and check data
			if (fem.Init() == false) return 1;
		}

		// solve the problem
		if (ops.brun)
		{
			bool bret = fem.Solve();
			return (int) (bret?0:1);
		}
		else
		{
			printf("\nData check complete\n\n");
		}
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : ParseCmdLine
//  Parses the command line and returns a CMDOPTIONS structure
//

bool ParseCmdLine(int nargs, char* argv[], CMDOPTIONS& ops)
{
	// set default options
	ops.brstrt = false;
	ops.boptim = false;
	ops.bdebug = false;
	ops.bdmp = false;
	ops.bplt = false;
	ops.blog = false;
	ops.brun = true;
	ops.bdiag = false;
	ops.bsplash = true;

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
			ops.brstrt = true;
			strcpy(ops.szfile, argv[++i]);
		}
		else if (strcmp(sz, "-d") == 0)
		{
			ops.bdiag = true;
			strcpy(ops.szfile, argv[++i]);
		}
		else if (strcmp(sz, "-p") == 0)
		{
			ops.bplt = true;
			strcpy(ops.szplt, argv[++i]);
		}
		else if (strcmp(sz, "-a") == 0)
		{
			ops.bdmp = true;
			strcpy(ops.szdmp, argv[++i]);
		}
		else if (strcmp(sz, "-o") == 0)
		{
			ops.blog = true;
			strcpy(ops.szlog, argv[++i]);
		}
		else if (strcmp(sz, "-i") == 0)
		{
			strcpy(ops.szfile, argv[++i]);
		}
		else if (strcmp(sz, "-s") == 0)
		{
			ops.boptim = true;
			strcpy(ops.szfile, argv[++i]);
		}
		else if (strcmp(sz, "-g") == 0)
		{
			ops.bdebug = true;
		}
		else if (strcmp(sz, "-c") == 0)
		{
			// don't run the problem.
			// just do a data check
			ops.brun = false;
		}
		else if (strcmp(sz, "-nosplash") == 0)
		{
			// don't show the welcome message
			ops.bsplash = false;
		}
		else if (strcmp(sz, "-cnf") == 0)
		{
			strcpy(ops.szcnf, argv[++i]);
		}
		else if (strcmp(sz, "-noconfig") == 0)
		{
			ops.szcnf[0] = 0;
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
		if (!ops.blog) sprintf(ops.szlog, "%s.log", szbase);
		if (!ops.bplt) sprintf(ops.szplt, "%s.plt", szbase);
		if (!ops.bdmp) sprintf(ops.szdmp, "%s.dmp", szbase);
	}
	else
	{
		if (!ops.blog) sprintf(ops.szlog, "%sf3log.txt", szpath);
		if (!ops.bplt) sprintf(ops.szplt, "%sf3plot"   , szpath);
		if (!ops.bdmp) sprintf(ops.szdmp, "%sf3dmp"    , szpath);
	}

	return true;
}

//-----------------------------------------------------------------------------
void init_framework(FEM& fem)
{
	FEBioCommand::SetFEM(&fem);
}

//-----------------------------------------------------------------------------
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
