/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEBioApp.h"
#include "console.h"
#include "CommandManager.h"
#include <FECore/log.h>
#include "console.h"
#include "breakpoint.h"
#include <FEBioLib/febio.h>
#include <FEBioLib/version.h>
#include "febio_cb.h"
#include "Interrupt.h"
#include "ping.h"

FEBioApp* FEBioApp::m_This = nullptr;

FEBioApp::FEBioApp()
{
	assert(m_This == nullptr);
	m_This = this;

	m_fem = nullptr;
}

FEBioApp* FEBioApp::GetInstance()
{
	assert(m_This);
	return m_This;
}

febio::CMDOPTIONS& FEBioApp::CommandOptions()
{
	return m_ops;
}

bool FEBioApp::Init(int argc, char* argv[])
{
	// Initialize kernel
	FECoreKernel::SetInstance(febio::GetFECoreKernel());

	// parse the command line
	if (ParseCmdLine(argc, argv) == false) return false;

	// say hello
	ConsoleStream s;
	if (m_ops.bsplash && (!m_ops.bsilent)) febio::Hello(s);

	// Initialize FEBio library
	febio::InitLibrary();

	// copy some flags to configuration
	m_config.SetOutputLevel(m_ops.bsilent ? 0 : 1);

	// read the configration file if specified
	if (m_ops.szcnf[0])
		if (febio::Configure(m_ops.szcnf, m_config) == false)
		{
			fprintf(stderr, "FATAL ERROR: An error occurred reading the configuration file.\n");
			return false;
		}

	// read command line plugin if specified
	if (m_ops.szimp[0] != 0)
	{
		febio::ImportPlugin(m_ops.szimp);
	}

	// ping repo server
	// Removed for the time being, pending further instruction
	// ping();

	return true;
}

//-----------------------------------------------------------------------------
bool FEBioApp::Configure(const char* szconfig)
{
	return febio::Configure(szconfig, m_config);
}

//-----------------------------------------------------------------------------
int FEBioApp::Run()
{
	// activate interruption handler
	Interruption I;

	// run FEBio either interactively or directly
	if (m_ops.binteractive)
		return prompt();
	else
		return RunModel();
}

//-----------------------------------------------------------------------------
void FEBioApp::Finish()
{
	febio::FinishLibrary();

	Console::GetHandle()->CleanUp();
}

//-----------------------------------------------------------------------------
// get the current model
FEBioModel* FEBioApp::GetCurrentModel()
{
	return m_fem;
}

//-----------------------------------------------------------------------------
// set the currently active model
void FEBioApp::SetCurrentModel(FEBioModel* fem)
{
	m_fem = fem;
}

//-----------------------------------------------------------------------------
// Run an FEBio input file. 
int FEBioApp::RunModel()
{
	// create the FEBioModel object
	FEBioModel fem;
	SetCurrentModel(&fem);

	// add console stream to log file
	if (m_ops.bsilent == false)
		fem.GetLogFile().SetLogStream(new ConsoleStream);
	else
		Console::GetHandle()->Deactivate();

	// register callbacks
	fem.AddCallback(update_console_cb, CB_MAJOR_ITERS | CB_INIT | CB_SOLVED | CB_STEP_ACTIVE, 0);
	fem.AddCallback(interrupt_cb, CB_ALWAYS, 0);
	fem.AddCallback(break_point_cb, CB_ALWAYS, 0);

	// set options that were passed on the command line
	fem.SetDebugLevel(m_ops.ndebug);
	fem.SetDumpLevel(m_ops.dumpLevel);
	fem.SetDumpStride(m_ops.dumpStride);

	// set the output filenames
	fem.SetLogFilename(m_ops.szlog);
	fem.SetPlotFilename(m_ops.szplt);
	fem.SetDumpFilename(m_ops.szdmp);

	// read the input file if specified
	int nret = 0;
	if (m_ops.szfile[0])
	{
		// read the input file
		if (fem.Input(m_ops.szfile) == false) nret = 1;
		else
		{
			// apply configuration overrides
			ApplyConfig(fem);
		}
	}

	// solve the model with the task and control file
	if (nret == 0)
	{
		bool bret = febio::SolveModel(fem, m_ops.sztask, m_ops.szctrl);

		nret = (bret ? 0 : 1);
	}

	// reset the current model pointer
	SetCurrentModel(nullptr);

	return nret;
}

//-----------------------------------------------------------------------------
// apply configuration changes to model
void FEBioApp::ApplyConfig(FEBioModel& fem)
{
	if (m_config.m_printParams != -1)
	{
		fem.SetPrintParametersFlag(m_config.m_printParams != 0);
	}
	fem.ShowWarningsAndErrors(m_config.m_bshowErrors);
}

//-----------------------------------------------------------------------------
//! Prints the FEBio prompt. If the user did not enter anything on the command
//! line when running FEBio then commands can be entered at the FEBio prompt.
//! This function returns the command arguments as a CMDOPTIONS structure.
int FEBioApp::prompt()
{
	// get a pointer to the console window
	Console* pShell = Console::GetHandle();

	// set the title
	pShell->SetTitle("FEBio4");

	// process commands
	ProcessCommands();

	return 0;
}

//-----------------------------------------------------------------------------
//!  Parses the command line and returns a CMDOPTIONS structure
//
bool FEBioApp::ParseCmdLine(int nargs, char* argv[])
{
	febio::CMDOPTIONS& ops = m_ops;

	// set default options
	ops.ndebug = 0;
	ops.bsplash = true;
	ops.bsilent = false;
	ops.binteractive = true;

	// these flags indicate whether the corresponding file name
	// was defined on the command line. Otherwise, a default name will be generated.
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

	// set initial configuration file name
	if (ops.szcnf[0] == 0)
	{
		char szpath[1024] = { 0 };
		febio::get_app_path(szpath, 1023);
		sprintf(ops.szcnf, "%sfebio.xml", szpath);
	}

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
		else if (strncmp(sz, "-dump_stride", 12) == 0)
		{
			if (sz[12] == '=')
			{
				ops.dumpStride = atoi(sz + 13);
				if (ops.dumpStride < 1)
				{
					fprintf(stderr, "FATAL ERROR: invalid dump stride.\n");
					return false;
				}
			}
			else
			{
				fprintf(stderr, "FATAL ERROR: missing '=' after -dump_stride.\n");
				return false;
			}
		}
		else if (strncmp(sz, "-dump", 5) == 0)
		{
			ops.dumpLevel = FE_DUMP_MAJOR_ITRS;
			if (sz[5] == '=') ops.dumpLevel = atoi(sz + 6);
			if ((ops.dumpLevel < 0) || (ops.dumpLevel > 3))
			{
				fprintf(stderr, "FATAL ERROR: invalid restart level.\n");
				return false;
			}

			if (i<nargs - 1)
			{
				char* szi = argv[i + 1];
				if (szi[0] != '-')
				{
					// assume this is the name of the dump file
					strcpy(ops.szdmp, argv[++i]);
					bdmp = true;
				}
			}
		}
		else if (strcmp(sz, "-o") == 0)
		{
			blog = true;
			strcpy(ops.szlog, argv[++i]);
		}
		else if (strcmp(sz, "-i") == 0)
		{
			++i;
			const char* szext = strrchr(argv[i], '.');
			if (szext == 0)
			{
				// we assume a default extension of .feb if none is provided
				sprintf(ops.szfile, "%s.feb", argv[i]);
			}
			else strcpy(ops.szfile, argv[i]);
			ops.binteractive = false;
		}
		else if (strcmp(sz, "-s") == 0)
		{
			if (ops.sztask[0] != 0) { fprintf(stderr, "-s is incompatible with other command line option.\n"); return false; }
			strcpy(ops.sztask, "optimize");
			strcpy(ops.szctrl, argv[++i]);
		}
		else if ((strcmp(sz, "-g") == 0) || (strcmp(sz, "-g1") == 0))
		{
			ops.ndebug = 1;
		}
		else if (strcmp(sz, "-g2") == 0)
		{
			ops.ndebug = 2;
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

			add_break_point(szbuf);
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

			char* szver = febio::getVersionString();

#ifdef _DEBUG
			fprintf(fp, "FEBio version  = %s (DEBUG)\n", szver);
#else
			fprintf(fp, "FEBio version  = %s\n", szver);
#endif
			if (fp != stdout) fclose(fp);
		}
		else if (strcmp(sz, "-norun") == 0)
		{
			brun = false;
		}

		else if (strcmp(sz, "-import") == 0)
		{
			if ((i < nargs - 1) && (argv[i+1][0] != '-'))
				strcpy(ops.szimp, argv[++i]);
			else
			{
				fprintf(stderr, "FATAL ERROR: insufficient number of arguments for -import.\n");
				return false;
			}
		}
		else if (sz[0] == '-')
		{
			fprintf(stderr, "FATAL ERROR: Invalid command line option.\n");
			return false;
		}
		else
		{
			// if no input file is given yet, we'll assume this is the input file
			if (ops.szfile[0] == 0)
			{
				const char* szext = strrchr(sz, '.');
				if (szext == 0)
				{
					// we assume a default extension of .feb if none is provided
					sprintf(ops.szfile, "%s.feb", sz);
				}
				else
				{
					strcpy(ops.szfile, sz);
				}
				ops.binteractive = false;
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

void FEBioApp::ProcessCommands()
{
	// get a pointer to the console window
	Console* pShell = Console::GetHandle();

	// get the command manager
	CommandManager* CM = CommandManager::GetInstance();

	// enter command loop
	int nargs;
	char* argv[32];
	while (1)
	{
		// get a command from the shell
		pShell->GetCommand(nargs, argv);
		if (nargs > 0)
		{
			// find the command that has this name
			Command* pcmd = CM->Find(argv[0]);
			if (pcmd)
			{
				int nret = pcmd->run(nargs, argv);
				if (nret == 1) break;
			}
			else
			{
				printf("Unknown command: %s\n", argv[0]);
			}
		}
        else if (nargs == 0) break;

		// make sure to clear the progress on the console
		FEBioModel* fem = GetCurrentModel();
		if ((fem == nullptr) || fem->IsSolved()) pShell->SetProgress(0);
	}
}
