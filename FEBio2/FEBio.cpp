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
// at the University of Utah. All rights reserved.
// Copyright (c) 2006 - 2013
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
#include "FEBioLib/validate.h"
#include "console.h"
#include "FECore/log.h"
#include "FEBioStdSolver.h"
#include "FECore/FECoreKernel.h"
#include "Interrupt.h"
#include "plugin.h"
#include "FEBioXML/XMLReader.h"

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
void Hello();
int prompt(CMDOPTIONS& ops);
int get_app_path (char *pname, size_t pathsize);
extern void InitFEBioLibrary();
bool Configure(FEBioModel& fem, const char *szfile);

//-----------------------------------------------------------------------------
// we use the console to log output 
class ConsoleStream : public LogStream
{
public:
	void print(const char* sz) { printf(sz); }
};

//-----------------------------------------------------------------------------
// callback to update window title
void update_console_cb(FEModel* pfem, void* pd)
{
	FEBioModel& fem = dynamic_cast<FEBioModel&>(*pfem);

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

	// print progress in title bar
	if (nsteps > 1)
		pShell->SetTitle("(step %d/%d: %.f%%) %s - %s", fem.m_nStep+1, nsteps, f, fem.GetFileTitle(), (bdebug?"FEBio (debug mode)": "FEBio"));
	else
		pShell->SetTitle("(%.f%%) %s - %s", f, fem.GetFileTitle(), (bdebug?"FEBio (debug mode)": "FEBio"));
}

//-----------------------------------------------------------------------------
// callback for ctrl+c interruptions
void interrupt_cb(FEModel* pfem, void* pd)
{
	Interruption itr;
	if (itr.m_bsig)
	{
		itr.m_bsig = false;
		itr.interrupt();
	}
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

	// load the license file
	LoadLicenseFile();

	// print welcome message
	if (ops.bsplash && (!ops.bsilent)) Hello();

	// if silent mode only output to file
	if (ops.bsilent) felog.SetMode(Logfile::FILE_ONLY);

	// initialize FEBio library
	InitFEBioLibrary();

	// if there are no arguments, print the FEBio prompt
	if (argc == 1)
	{
		 if (prompt(ops) == 0) return 0;
	}

	// create the one and only FEBioModel object
	FEBioModel fem;

	// register callbacks
	fem.AddCallback(update_console_cb, CB_MAJOR_ITERS, 0);
	fem.AddCallback(interrupt_cb     , CB_MINOR_ITERS, 0);

	// intialize the framework
	FEBioCommand::SetFEM(&fem);

	// read the configration file if specified
	if (ops.szcnf[0])
	{
		if (Configure(fem, ops.szcnf) == false) return 1;
	}

	// set options that were passed on the command line
	fem.SetDebugFlag(ops.bdebug);

	// set the output filenames
	fem.SetLogFilename (ops.szlog);
	fem.SetPlotFilename(ops.szplt);
	fem.SetDumpFilename(ops.szdmp);

	// read the input file if specified
	if (ops.szfile[0])
	{
		// store the input file name
		fem.SetInputFilename(ops.szfile);

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

	// run the FEBio task (and pass the optional control file)
	bool bret = ptask->Run(ops.szctrl);

	// Don't forget to cleanup the plugins
	FEBioPluginManager* pPM = FEBioPluginManager::GetInstance();
	pPM->DeleteThis();

	// return the error code of the run
	return (bret?0:1);
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

	ops.szfile[0] = 0;
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
		else if (strcmp(sz, "-cnf") == 0)
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
		else if (strcmp(sz, "-import") == 0)
		{
			char* szfile = argv[++i];
			FEBioPluginManager* pPM = FEBioPluginManager::GetInstance();
			if (pPM->LoadPlugin(szfile) == false)
			{
				fprintf(stderr, "Failed loading plugin %s\n\n", szfile);
			}
			else
			{
				fprintf(stderr, "Success loading plugin %s\n\n", szfile);
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
					fprintf(stderr, "FATAL ERROR: Invalid command line option\n\n");
					return false;
				}
			}
			else
			{
				fprintf(stderr, "FATAL ERROR: Invalid command line option\n\n");
				return false;
			}
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

		if (!blog) sprintf(ops.szlog, "%s.log", szbase);
		if (!bplt) sprintf(ops.szplt, "%s.xplt", szbase);
		if (!bdmp) sprintf(ops.szdmp, "%s.dmp", szbase);
	}

	return brun;
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

//-----------------------------------------------------------------------------
//! Reads the FEBio configuration file. This file contains some default settings.

bool Configure(FEBioModel& fem, const char *szfile)
{
	// open the configuration file
	XMLReader xml;
	if (xml.Open(szfile) == false)
	{
		fprintf(stderr, "FATAL ERROR: Failed reading FEBio configuration file %s.\n", szfile);
		return false;
	}

	// loop over all child tags
	try
	{
		// Find the root element
		XMLTag tag;
		if (xml.FindTag("febio_config", tag) == false) return false;

		if (strcmp(tag.m_att[0].m_szatv, "1.0") == 0)
		{
			if (!tag.isleaf())
			{
				// Read version 1.0
				++tag;
				do
				{
					if (tag == "linear_solver")
					{
						const char* szt = tag.AttributeValue("type");
						if      (strcmp(szt, "skyline"           ) == 0) fem.m_nsolver = SKYLINE_SOLVER;
						else if (strcmp(szt, "psldlt"            ) == 0) fem.m_nsolver = PSLDLT_SOLVER;
						else if (strcmp(szt, "superlu"           ) == 0) fem.m_nsolver = SUPERLU_SOLVER;
						else if (strcmp(szt, "superlu_mt"        ) == 0) fem.m_nsolver = SUPERLU_MT_SOLVER;
						else if (strcmp(szt, "pardiso"           ) == 0) fem.m_nsolver = PARDISO_SOLVER;
						else if (strcmp(szt, "rcicg"             ) == 0) fem.m_nsolver = RCICG_SOLVER;
						else if (strcmp(szt, "wsmp"              ) == 0) fem.m_nsolver = WSMP_SOLVER;
					}
					else if (tag == "import")
					{
						const char* szfile = tag.szvalue();
						FEBioPluginManager* pPM = FEBioPluginManager::GetInstance();
						if (pPM->LoadPlugin(szfile) == false) throw XMLReader::InvalidValue(tag);
						printf("Plugin \"%s\" loaded successfully\n", szfile);
					}
/*					else if (tag == "import_folder")
					{
						const char* szfile = tag.szvalue();
						if (LoadPluginFolder(szfile) == false) throw XMLReader::InvalidTag(tag);
					}
*/					else if (tag == "omp_num_threads")
					{
						int n;
						tag.value(n);
						  //						omp_set_num_threads(n);
					}
					else throw XMLReader::InvalidTag(tag);

					// go to the next tag
					++tag;
				}
				while (!tag.isend());
			}
		}
		else
		{
			felog.printbox("FATAL ERROR", "Invalid version for FEBio configuration file.");
			return false;
		}
	}
	catch (XMLReader::XMLSyntaxError)
	{
		felog.printf("FATAL ERROR: Syntax error (line %d)\n", xml.GetCurrentLine());
		return false;
	}
	catch (XMLReader::InvalidTag e)
	{
		felog.printf("FATAL ERROR: unrecognized tag \"%s\" (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::InvalidAttributeValue e)
	{
		const char* szt = e.tag.m_sztag;
		const char* sza = e.szatt;
		const char* szv = e.szval;
		int l = e.tag.m_nstart_line;
		felog.printf("FATAL ERROR: unrecognized value \"%s\" for attribute \"%s.%s\" (line %d)\n", szv, szt, sza, l);
		return false;
	}
	catch (XMLReader::InvalidValue e)
	{
		felog.printf("FATAL ERROR: the value for tag \"%s\" is invalid (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::MissingAttribute e)
	{
		felog.printf("FATAL ERROR: Missing attribute \"%s\" of tag \"%s\" (line %d)\n", e.szatt, e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::UnmatchedEndTag e)
	{
		const char* sz = e.tag.m_szroot[e.tag.m_nlevel];
		felog.printf("FATAL ERROR: Unmatched end tag for \"%s\" (line %d)\n", sz, e.tag.m_nstart_line);
		return false;
	}
	catch (...)
	{
		felog.printf("FATAL ERROR: unrecoverable error (line %d)\n", xml.GetCurrentLine());
		return false;
	}

	xml.Close();

	return true;
}
