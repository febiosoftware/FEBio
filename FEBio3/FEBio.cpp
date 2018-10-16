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
#include "console.h"
#include "FECore/log.h"
#include "FEBioStdSolver.h"
#include "FECore/FECoreKernel.h"
#include "FECore/FEAnalysis.h"
#include "Interrupt.h"
#include "FEBioXML/XMLReader.h"
#include <FEBioLib/febio.h>
#include <FEBioLib/plugin.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef WIN32
extern "C" void __cdecl omp_set_num_threads(int);
#else
extern "C" void omp_set_num_threads(int);
#endif

//-----------------------------------------------------------------------------
// TODO: On Windows the GetCurrentTime macro gets in here via plugin.h. 
// I need to look into how to prevent this
#ifdef GetCurrentTime
#undef GetCurrentTime
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
	bool	binteractive;		//!< start FEBio interactively

	char	szfile[MAXFILE];	//!< model input file name
	char	szlog [MAXFILE];	//!< log file name
	char	szplt [MAXFILE];	//!< plot file name
	char	szdmp [MAXFILE];	//!< dump file name
	char	szcnf [MAXFILE];	//!< configuration file
	char	sztask[MAXFILE];	//!< task name
	char	szctrl[MAXFILE];	//!< control file for tasks
	char	szimp[MAXFILE];		//!< import file
};

//-----------------------------------------------------------------------------
// forward declarations
//
bool ParseCmdLine(int argc, char* argv[], CMDOPTIONS& ops);
int Run(CMDOPTIONS& ops);
int prompt(CMDOPTIONS& ops);

//-----------------------------------------------------------------------------
// we use the console to log output 
class ConsoleStream : public LogStream
{
public:
	void print(const char* sz) { printf("%s",sz); }
};

//-----------------------------------------------------------------------------
// callback to update window title
bool update_console_cb(FEModel* pfem, unsigned int nwhen, void* pd)
{
	FEBioModel& fem = static_cast<FEBioModel&>(*pfem);

	// get the number of steps
	int nsteps = fem.Steps();

	// calculate progress
	double starttime = fem.GetStartTime();
	double endtime = fem.GetCurrentStep()->m_tend;
	double f = 0.0;
	double ftime = fem.GetCurrentTime();
	if (endtime > 0.0) f = 100.f*(ftime - starttime) / (endtime - starttime);

	// check debug flag
	bool bdebug = fem.GetDebugFlag();

	// obtain a pointer to the console object. We'll use this to
	// set the title of the console window.
	Console* pShell = Console::GetHandle();

	char szvers[32] = {0};
#ifdef _DEBUG
	sprintf(szvers, "FEBio (DEBUG) %d.%d.%d", VERSION, SUBVERSION, SUBSUBVERSION);
#else
	sprintf(szvers, "FEBio %d.%d.%d", VERSION, SUBVERSION, SUBSUBVERSION);
#endif

	// print progress in title bar
	const char* szfile = fem.GetFileTitle();
	if (szfile == 0) szfile = "";

	if (nsteps > 1)
		pShell->SetTitle("(step %d/%d: %.f%%) %s - %s %s", fem.GetCurrentStepIndex() + 1, nsteps, f, szfile, szvers, (bdebug?"(debug mode)": ""));
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
	double t = pfem->GetTime().currentTime;

	Interruption itr;
	int nbp = (int)break_points.size();
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
#ifdef USE_MPI
	MPI_Init(&argc, &argv);
#endif

	// divert the log output to the console
	felog.SetLogStream(new ConsoleStream);

	// parse the command line
	CMDOPTIONS ops;
	if (ParseCmdLine(argc, argv, ops) == false) return 0;

	// say hello
	if (ops.bsplash && (!ops.bsilent)) febio::Hello();

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

	// read command line plugin if specified
	if (ops.szimp[0] != 0)
	{
		febio::ImportPlugin(ops.szimp);
	}

	// start FEBio
	int nret = 0;
	if (ops.binteractive)	 
		nret = (prompt(ops));
	else
		nret = Run(ops);

	// Don't forget to cleanup
	febio::FinishLibrary();

#ifdef USE_MPI
	MPI_Finalize();
#endif

	return nret;
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
	ops.binteractive = true;

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
	ops.szimp[0] = 0;

	// set the location of the configuration file
	char szpath[1024] = {0};
	febio::get_app_path (szpath, 1023);

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
			ops.binteractive = false;
		}
		else if (strcmp(sz, "-d") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-d is incompatible with other command line option.\n"); return false; }
			strcpy(ops.sztask, "diagnose");
			strcpy(ops.szctrl, argv[++i]);
			ops.binteractive = false;
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
			ops.binteractive = false;
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
#ifdef _DEBUG
			fprintf(fp, "FEBio version  = %d.%d.%d (DEBUG)\n", VERSION, SUBVERSION, SUBSUBVERSION);
#else
			fprintf(fp, "FEBio version  = %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
#endif
			fprintf(fp, "FECore version = %s\n", FECore::get_version_string());
			if (fp != stdout) fclose(fp);
		}
		else if (strcmp(sz, "-norun") == 0)
		{
			brun = false;
		}

		else if (strcmp(sz, "-import") == 0)
		{
			strcpy(ops.szimp, argv[++i]);
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
					ops.binteractive = false;
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
		char* ch = strrchr(szbase, '.');
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
		char* ch = strrchr(szbase, '.');
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
		felog.SetMode(Logfile::LOG_FILE);
		Console::GetHandle()->Deactivate();
	}

	// create the one and only FEBioModel object
	FEBioModel fem;

	// register callbacks
	fem.AddCallback(update_console_cb, CB_MAJOR_ITERS | CB_INIT, 0);
	fem.AddCallback(interrupt_cb     , CB_ALWAYS, 0);
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
	FECoreTask* ptask = fecore_new<FECoreTask>(ops.sztask, &fem);
	if (ptask == 0)
	{
		fprintf(stderr, "Don't know how to do task: %s\n", ops.sztask);
		return 1;
	}

	// initialize the task
	if (ptask->Init(ops.szctrl) == false)
	{
		fprintf(stderr, "Failed initializing the task: %s\n", ops.sztask);
		return 1;
	}

	// run the task
	bool bret = ptask->Run();

	return (bret?0:1);
}

//-----------------------------------------------------------------------------
void cmd_help()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "config [file]   ...............   (re-)load a FEBio configuration file\n");
	fprintf(stderr, "help  .........................   print this info\n");
	fprintf(stderr, "load [file]  ..................   load a plugin\n");
	fprintf(stderr, "plugins  ......................   list the plugins that are loaded\n");
	fprintf(stderr, "quit  .........................   exits FEBio\n");
	fprintf(stderr, "run [-i,-s] <file> [OPTIONS] ..   run an FEBio input file\n");
	fprintf(stderr, "unload [file|%%n] .............   unload a plugin\n");
	fprintf(stderr, "version  ......................   print version information\n");
}

//-----------------------------------------------------------------------------
void cmd_run(int nargs, char* argv[], CMDOPTIONS& ops)
{
	ParseCmdLine(nargs, argv, ops);

	// run the FEBio3 on the ops
	Run(ops);

	// reset the title after computation.
	Console* pShell = Console::GetHandle();
	pShell->SetTitle("FEBio3");
}

//-----------------------------------------------------------------------------
void cmd_load(int nargs, char* argv[])
{
	if (nargs < 2) fprintf(stderr, "missing file name\n");
	else febio::ImportPlugin(argv[1]);
}

//-----------------------------------------------------------------------------
void cmd_unload(int nargs, char* argv[])
{
	FEBioPluginManager* PM = FEBioPluginManager::GetInstance();
	if (PM == 0) return;

	if (nargs == 1)
	{
		// unload all plugins
		while (PM->Plugins() > 0)
		{
			const FEBioPlugin& pl = PM->GetPlugin(0);
			string sname = pl.GetName();
			bool b = PM->UnloadPlugin(0);
			if (b) fprintf(stdout, "Success unloading %s\n", sname.c_str());
			else fprintf(stdout, "Failed unloading %s\n", sname.c_str());
		}
	}
	else if (nargs == 2)
	{
		const char* c = argv[1];
		if (c[0] == '%')
		{
			int n = atoi(c+1);
			if ((n > 0) && (n <= PM->Plugins()))
			{
				const FEBioPlugin& pl = PM->GetPlugin(n - 1);
				string sname = pl.GetName();
				bool b = PM->UnloadPlugin(n - 1);
				if (b) fprintf(stdout, "Success unloading %s\n", sname.c_str());
				else fprintf(stdout, "Failed unloading %s\n", sname.c_str());
			}
			else fprintf(stderr, "Invalid plugin index\n");
		}
		else
		{
			bool b = PM->UnloadPlugin(argv[1]);
			if (b) fprintf(stdout, "Success unloading %s\n", argv[1]);
			else fprintf(stdout, "Failed unloading %s\n", argv[1]);
		}
	}
	else fprintf(stderr, "syntax error\n");
}

//-----------------------------------------------------------------------------
void cmd_version()
{
#ifdef _DEBUG
	fprintf(stderr, "\nFEBio version %d.%d.%d (DEBUG)\n", VERSION, SUBVERSION, SUBSUBVERSION);
#else
	fprintf(stderr, "\nFEBio version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
#endif
	fprintf(stderr, "SDK Version %d.%d\n", FE_SDK_MAJOR_VERSION, FE_SDK_SUB_VERSION);
	fprintf(stderr, "FECore version %s\n", FECore::get_version_string());
	fprintf(stderr, "compiled on " __DATE__ "\n\n");
}

//-----------------------------------------------------------------------------
void cmd_config(int nargs, char* argv[], CMDOPTIONS& ops)
{
	if (nargs == 1)
	{
		febio::Configure(ops.szcnf);
	}
	else if (nargs == 2)
	{
		char szpath[1024] = { 0 };
		febio::get_app_path(szpath, 1023);
		sprintf(ops.szcnf, "%s%s", szpath, argv[1]);
		febio::Configure(ops.szcnf);
	}
	else
	{
		printf("Invalid number of command arguments.\n");
	}
}

//-----------------------------------------------------------------------------
void cmd_plugins()
{
	FEBioPluginManager* PM = FEBioPluginManager::GetInstance();
	if (PM == 0) return;

	int NP = PM->Plugins();
	for (int i=0; i<NP; ++i)
	{
		const FEBioPlugin& pl = PM->GetPlugin(i);
		fprintf(stdout, "%%%d: %s\n", i+1, pl.GetName());
	}
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
	pShell->SetTitle("FEBio3");
	int nargs;
	char* argv[32];

	while (1)
	{
		// get a command from the shell
		pShell->GetCommand(nargs, argv);
		if (nargs > 0)
		{
			if (strcmp(argv[0], "quit") == 0) return 0;
			else if (strcmp(argv[0], "help"   ) == 0) cmd_help();
			else if (strcmp(argv[0], "run"    ) == 0) cmd_run(nargs, argv, ops);
			else if (strcmp(argv[0], "load"   ) == 0) cmd_load(nargs, argv);
			else if (strcmp(argv[0], "unload" ) == 0) cmd_unload(nargs, argv);
			else if (strcmp(argv[0], "version") == 0) cmd_version();
			else if (strcmp(argv[0], "config" ) == 0) cmd_config(nargs, argv, ops);
			else if (strcmp(argv[0], "plugins") == 0) cmd_plugins();
			else
			{
				printf("Unknown command: %s\n", argv[0]);
				fprintf(stderr, "Type help for an overview of commands.\n");
			}
		}
	}
	return 0;
}
