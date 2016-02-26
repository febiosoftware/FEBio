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
// FEBio is a finite element solver that is specifically designed for three
// dimensional biomechanical applications. It solves the nonlinear finite
// element equations using a quasi-Newton method called the BFGS-method. It
// also offers several biologically relevant constitutive models.
//
// This software is developed at the Musculoskeletal Research Laboratories
// at the University of Utah. FEBio is a registered trademark. All rights reserved. 
// Copyright (c) 2006 - 2015
//
// The subversion (svn) revision number of this code can be found in the file
// FEBio/svnrev.h
//
// Main developers:
//  - Steve Maas
//  - Jeff Weiss
//  - Dave Rawlins
//  - Gerard Ateshian
//
// Contributors:
//  - Alexander Veress
//
///////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEBioLib/FEBioModel.h"
#include "FEBioLib/version.h"
#include "FEBioCommand.h"
#include "FECore/FECore.h"
//#include "FEBioLib/validate.h" // For the KeyGen library
#include "console.h"
#include "FECore/log.h"
#include "FEBioStdSolver.h"
#include "FECore/FECoreKernel.h"
#include "Interrupt.h"
#include "FEBioXML/XMLReader.h"
#include <FEBioLib/febio.h>
#include <FEBioLib/plugin.h>

#ifdef WIN32
extern "C" void __cdecl omp_set_num_threads(int);
#else
//extern "C" void omp_set_num_threads(int);
#endif

//-----------------------------------------------------------------------------
//!  Command line options

//! This structures stores the command line options that were input by the user
struct CMDOPTIONS
{
	enum {MAXFILE=512};

	bool		bdebug;			//!< debug flag

	bool	bsplash;			//!< show splash screen or not
	bool	bsilent;			//!< run FEBio in silent mode (no output to screen)

	char	szfile[MAXFILE];	//!< model input file name
	char	szlog [MAXFILE];	//!< log file name
	char	szplt [MAXFILE];	//!< plot file name
	char	szdmp [MAXFILE];	//!< dump file name
	char	szcnf [MAXFILE];	//!< configuration file
	char	sztask[MAXFILE];	//!< task name
	char	szctrl[MAXFILE];	//!< control file for tasks
};

//-----------------------------------------------------------------------------
// forward declarations
//
bool ParseCmdLine(int argc, char* argv[], CMDOPTIONS& ops);
int Hello();
int Run(CMDOPTIONS& ops);
int prompt(CMDOPTIONS& ops);
int get_app_path (char *pname, size_t pathsize);

//-----------------------------------------------------------------------------
// we use the console to log output 
class ConsoleStream : public LogStream
{
public:
	void print(const char* sz) { printf(sz); }
};

//-----------------------------------------------------------------------------
// callback to update window title
bool update_console_cb(FEModel* pfem, unsigned int nwhen, void* pd)
{
	FEBioModel& fem = static_cast<FEBioModel&>(*pfem);

	// get the number of steps
	int nsteps = fem.Steps();

	// calculate progress
	double starttime = fem.m_ftime0;
	double endtime = fem.GetCurrentStep()->m_tend;
	double f = 0.0;
	if (endtime > 0.0) f = 100.f*(fem.m_ftime - starttime) / (endtime - starttime);

	// check debug flag
	bool bdebug = fem.GetDebugFlag();

	// obtain a pointer to the console object. We'll use this to
	// set the title of the console window.
	Console* pShell = Console::GetHandle();

	char szvers[32] = {0};
#ifdef _DEBUG
	sprintf(szvers, "FEBio (DEBUG BUILD) %d.%d.%d.%d", VERSION, SUBVERSION, SUBSUBVERSION, SVNREVISION);
#else
	sprintf(szvers, "FEBio %d.%d.%d.%d", VERSION, SUBVERSION, SUBSUBVERSION, SVNREVISION);
#endif

	// print progress in title bar
	const char* szfile = fem.GetFileTitle();
	if (szfile == 0) szfile = "";

	if (nsteps > 1)
		pShell->SetTitle("(step %d/%d: %.f%%) %s - %s %s", fem.m_nStep+1, nsteps, f, szfile, szvers, (bdebug?"(debug mode)": ""));
	else
		pShell->SetTitle("(%.f%%) %s - %s %s", f, szfile, szvers, (bdebug?"(debug mode)": ""));

	return true;
}

//-----------------------------------------------------------------------------
vector<pair<double,int> > break_points;

void clear_break_points()
{
	break_points.clear();
}

void add_break_point(double t)
{
	pair<double,int> bp;
	bp.first = t;
	bp.second = 1;
	break_points.push_back(bp);
}

//-----------------------------------------------------------------------------
// break points cb
bool break_point_cb(FEModel* pfem, unsigned int nwhen, void* pd)
{
	const double eps = 1e-12;

	// see if a break-point has been reached
	double t = pfem->GetTime().t;

	Interruption itr;
	int nbp = break_points.size();
	for (int i=0; i<nbp; ++i)
	{
		pair<double, int>& bpi = break_points[i];
		if (bpi.second && (bpi.first <= t + eps))
		{
			itr.interrupt();
			bpi.second = 0;
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
// callback for ctrl+c interruptions
bool interrupt_cb(FEModel* pfem, unsigned int nwhen, void* pd)
{
	Interruption itr;
	if (itr.m_bsig)
	{
		itr.m_bsig = false;
		itr.interrupt();
	}
	return true;
}

//-----------------------------------------------------------------------------
// The starting point of the application
//
int main(int argc, char* argv[])
{
	// divert the log output to the console
	felog.SetLogStream(new ConsoleStream);

	// parse the command line
	CMDOPTIONS ops;
	if (ParseCmdLine(argc, argv, ops) == false) return 0;

	// load the license file (for KeyGen)
	//LoadLicenseFile();

	// print welcome message
#ifdef NALPLIB
	if (Hello()) return 1;
#else
	if (ops.bsplash && (!ops.bsilent)) Hello();
#endif
	// Initialize FEBio library
	febio::InitLibrary();

	// set default linear solver
	// (Set this before the configuration is read in because
	//  the configuration can change the default linear solver.)
#ifdef PARDISO
	FECoreKernel::SetDefaultSolver(PARDISO_SOLVER);
#else
	FECoreKernel::SetDefaultSolver(SKYLINE_SOLVER);
#endif

	// read the configration file if specified
	if (ops.szcnf[0])
		if (febio::Configure(ops.szcnf) == false) return 1;

	// if there are no arguments, print the FEBio prompt
	if (argc == 1)	 return (prompt(ops));
	else			 return Run(ops);

	// Don't forget to cleanup
	febio::FinishLibrary();
}

//-----------------------------------------------------------------------------
//!  Parses the command line and returns a CMDOPTIONS structure
//
bool ParseCmdLine(int nargs, char* argv[], CMDOPTIONS& ops)
{
	// set default options
	ops.bdebug = false;
	ops.bsplash = true;
	ops.bsilent = false;

	bool blog = false;
	bool bplt = false;
	bool bdmp = false;
	bool brun = true;

	// initialize file names
	ops.szfile[0] = 0;
	ops.szplt[0] = 0;
	ops.szlog[0] = 0;
	ops.szdmp[0] = 0;
	ops.sztask[0] = 0;
	ops.szctrl[0] = 0;

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
			if (ops.sztask[0] != 0) { fprintf(stderr, "-r is incompatible with other command line option.\n"); return false; }
			strcpy(ops.sztask, "restart");
			strcpy(ops.szctrl, argv[++i]);
		}
		else if (strcmp(sz, "-d") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-d is incompatible with other command line option.\n"); return false; }
			strcpy(ops.sztask, "diagnose");
			strcpy(ops.szctrl, argv[++i]);
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
			strcpy(ops.szfile, argv[++i]);
		}
		else if (strcmp(sz, "-s") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-s is incompatible with other command line option.\n"); return false; }
			strcpy(ops.sztask, "optimize");
			strcpy(ops.szctrl, argv[++i]);
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
		else if (strcmp(sz, "-cnf") == 0)	// obsolete: use -config instead
		{
			strcpy(ops.szcnf, argv[++i]);
		}
		else if (strcmp(sz, "-config") == 0)
		{
			strcpy(ops.szcnf, argv[++i]);
		}
		else if (strcmp(sz, "-noconfig") == 0)
		{
			ops.szcnf[0] = 0;
		}
		else if (strncmp(sz, "-task", 5) == 0)
		{
			if (sz[5] != '=') { fprintf(stderr, "command line error when parsing task\n"); return false; }
			strcpy(ops.sztask, sz+6);

			if (i<nargs-1)
			{
				char* szi = argv[i+1];
				if (szi[0] != '-')
				{
					// assume this is a control file for the specified task
					strcpy(ops.szctrl, argv[++i]);
				}
			}
		}
		else if (strcmp(sz, "-break") == 0)
		{
			char szbuf[32]={0};
			strcpy(szbuf, argv[++i]);
			double f = atof(szbuf);
			add_break_point(f);
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
			if (SVNREVISION) fprintf(fp, "SVN revision: %d\n", SVNREVISION);
			fprintf(fp, "FECore version = %s\n", FECore::get_version_string());
			fprintf(fp, "SVN revision   = %d\n", SVNREVISION);
			if (fp != stdout) fclose(fp);
		}
		else if (strcmp(sz, "-norun") == 0)
		{
			brun = false;
		}

		else if (strcmp(sz, "-import") == 0)
		{
			char* szfile = argv[++i];
			FEBioPluginManager* pPM = FEBioPluginManager::GetInstance();
			int nerr = pPM->LoadPlugin(szfile);
			switch (nerr)
			{
			case 0: fprintf(stderr, "Success loading plugin %s\n", szfile); break;
			case 1: fprintf(stderr, "Failed loading plugin %s\n Reason: Failed to load the file.\n\n", szfile); break;
			case 2: fprintf(stderr, "Failed loading plugin %s\n Reason: Required plugin function PluginNumClasses not found.\n\n", szfile); break;
			case 3: fprintf(stderr, "Failed loading plugin %s\n Reason: Required plugin function PluginGetFactory not found.\n\n", szfile); break;
			case 4: fprintf(stderr, "Failed loading plugin %s\n Reason: Invalid number of classes returned by PluginNumClasses.\n\n", szfile); break;
			default:
				fprintf(stderr, "Failed loading plugin %s\n Reason: unspecified.\n\n", szfile); break;
			}
		}
		else
		{
			// we allow FEBio to run without a -i option
			// so that we can run an .feb file by right-clicking on it in windows
			if (nargs == 2)
			{
				char* c = strrchr(argv[1], '.');
				if (c && (strcmp(c, ".feb")==0))
				{
					strcpy(ops.szfile, argv[1]);
				}
				else
				{
					fprintf(stderr, "FATAL ERROR: Invalid command line option\n");
					return false;
				}
			}
			else
			{
				fprintf(stderr, "FATAL ERROR: Invalid command line option\n");
				return false;
			}
		}
	}

	// do some sanity checks
	if (strcmp(ops.sztask, "optimize") == 0)
	{
		// make sure we have an input file
		if (ops.szfile[0]==0)
		{
			fprintf(stderr, "FATAL ERROR: no model input file was defined (use -i to define the model input file)\n\n");
			return false;
		}
	}

	// if no task is defined, we assume a std solve is wanted
	if (ops.sztask[0] == 0) strcpy(ops.sztask, "solve");

	// derive the other filenames
	if (ops.szfile[0])
	{
		char szbase[256]; strcpy(szbase, ops.szfile);
		ch = strrchr(szbase, '.');
		if (ch) *ch = 0;

		char szlogbase[256];
		if (ops.szctrl[0])
		{
			strcpy(szlogbase, ops.szctrl);
			ch = strrchr(szlogbase, '.');
			if (ch) *ch = 0;
		}
		else strcpy(szlogbase, szbase);

		if (!blog) sprintf(ops.szlog, "%s.log", szlogbase);
		if (!bplt) sprintf(ops.szplt, "%s.xplt", szbase);
		if (!bdmp) sprintf(ops.szdmp, "%s.dmp", szbase);
	}
	else if (ops.szctrl[0])
	{
		char szbase[256]; strcpy(szbase, ops.szfile);
		strcpy(szbase, ops.szctrl);
		ch = strrchr(szbase, '.');
		if (ch) *ch = 0;

		if (!blog) sprintf(ops.szlog, "%s.log", szbase);
		if (!bplt) sprintf(ops.szplt, "%s.xplt", szbase);
		if (!bdmp) sprintf(ops.szdmp, "%s.dmp", szbase);
	}


	return brun;
}

//-----------------------------------------------------------------------------
// Run an FEBio input file. 
int Run(CMDOPTIONS& ops)
{
	// if silent mode only output to file
	if (ops.bsilent)
	{
		felog.SetMode(Logfile::FILE_ONLY);
		Console::GetHandle()->Deactivate();
	}

	// create the one and only FEBioModel object
	FEBioModel fem;

	// register callbacks
	fem.AddCallback(update_console_cb, CB_MAJOR_ITERS | CB_INIT, 0);
	fem.AddCallback(interrupt_cb     , CB_MINOR_ITERS, 0);
	fem.AddCallback(break_point_cb   , CB_MAJOR_ITERS, 0);

	// intialize the framework
	FEBioCommand::SetFEM(&fem);

	// set options that were passed on the command line
	fem.SetDebugFlag(ops.bdebug);

	// set the output filenames
	fem.SetLogFilename (ops.szlog);
	fem.SetPlotFilename(ops.szplt);  
	fem.SetDumpFilename(ops.szdmp);
		
	// read the input file if specified
	if (ops.szfile[0])
	{
		// read the input file
		if (fem.Input(ops.szfile) == false) return 1;
	}

	// find a task
	FECoreTask* ptask = fecore_new<FECoreTask>(FETASK_ID, ops.sztask, &fem);
	if (ptask == 0)
	{
		fprintf(stderr, "Don't know how to do task: %s\n", ops.sztask);
		return 1;
	}

	// initialize the task
	if (ptask->Init(ops.szctrl) == false)
	{
		fprintf(stderr, "Failed initializing the task: %s", ops.sztask);
		return 1;
	}

	// run the task
	bool bret = ptask->Run();

	return (bret?0:1);
}

//-----------------------------------------------------------------------------
//! Prints the FEBio prompt. If the user did not enter anything on the command
//! line when running FEBio then commands can be entered at the FEBio prompt.
//! This function returns the command arguments as a CMDOPTIONS structure.
int prompt(CMDOPTIONS& ops)
{
	// get a pointer to the console window
	Console* pShell = Console::GetHandle();

	// set the title
	pShell->SetTitle("FEBio2");
	int nargs;
	char* argv[32];

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
				fprintf(stderr, "import - load a plugin\n");
				fprintf(stderr, "quit - exits the application\n");
				fprintf(stderr, "run [-i,-s] <file> [OPTIONS] - run an FEBio input file\n");
				fprintf(stderr, "version - print version information\n");
			}
			else if (strcmp(argv[0], "run") == 0)
			{
				ParseCmdLine(nargs, argv, ops);

				// run the FEBio2 on the ops
				Run(ops);
					
				// reset the title after computation.
				pShell->SetTitle("FEBio2");
			}
			else if (strcmp(argv[0], "import") == 0)
			{
				if (nargs < 2) fprintf(stderr, "missing file name\n");
				else febio::ImportPlugin(argv[1]);
			}
			else if (strcmp(argv[0], "version") == 0)
			{
#ifdef _WIN64
				fprintf(stderr, "\nFEBio version %d.%d.%d (x64)\n", VERSION, SUBVERSION, SUBSUBVERSION);
#else
				fprintf(stderr, "\nFEBio version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
#endif
				fprintf(stderr, "SVN revision: %d\n", SVNREVISION);
				fprintf(stderr, "compiled on " __DATE__ "\n");
				fprintf(stderr, "using FECore version %s\n\n", FECore::get_version_string());
			}
			else
			{
				printf("Unknown command: %s\n", argv[0]);
				fprintf(stderr, "Type help for an overview of commands.\n");
			}
		}
	}
	return 0;
}
