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
//! This software is developed at the Scientific Computing and Imaging institute
//! at the University of Utah. All rights reserved.
//! Copyright (c) 2006 - 2008
//! 
//! Main developers:	
//!  - Steve Maas
//!  - Jeff Weiss
//!
//! Contributors:
//!  - Alexander Veress
//!  - Gerard Ateshian
//
///////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FESolver.h"
#include "FEBioCommand.h"
#include "FECore/FECore.h"


#define MAXFILE 256

//-----------------------------------------------------------------------------
//!  Command line options

//! This structures stores the command line options that where input by the user

struct CMDOPTIONS
{
	char	szfile[MAXFILE];	//!< input file name

	bool	blog;
	char	szlog[MAXFILE];		//!< log file name

	bool	bplt;
	char	szplt[MAXFILE];		//!< plot file name

	bool	bdmp;				//!< dump flag
	char	szdmp[MAXFILE];		//!< dump file name

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

bool optimize(FEM& fem, const char* szfile);

bool diagnose(FEM& fem, const char* szfile);

void init_framework(FEM& fem);

//-----------------------------------------------------------------------------
// The starting point of the application

int main(int argc, char* argv[])
{
	// parse the command line
	CMDOPTIONS ops;
	if (ParseCmdLine(argc, argv, ops) == false) return 0;

	// print welcome message
	if (ops.bsplash) Hello(stdout);

	// create the one and only FEM object
	FEM fem;

	// intialize the framework
	init_framework(fem);


	// set the filenames
	fem.SetLogFilename (ops.szlog);
	fem.SetPlotFilename(ops.szplt);
	fem.SetDumpFilename(ops.szdmp);

	// set options that were passed on the command line
	fem.SetDebugFlag(ops.bdebug);

	// set the default plot and print levels
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

	// check for parameter optimization
	if (ops.boptim) optimize(fem, ops.szfile);
	else if (ops.bdiag) diagnose(fem, ops.szfile);
	else
	{
		// input data
		if (ops.brstrt)
		{
			// do a restart
			if (fem.Restart(ops.szfile) == false) return 0;
		}
		else
		{
			// read input data
			if (fem.Input(ops.szfile) == false) return 0;

			// initialize and check data 
			if (fem.Init() == false) return 0;
		}

		// solve the problem
		if (ops.brun) 
		{
			fem.Solve();
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

	// if there are no arguments, ask for an input file
	if (nargs == 1)
	{
		printf("Enter input file: ");
		scanf("%s", ops.szfile);
	}
	else
	{
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
			else
			{
				fprintf(stderr, "FATAL ERROR: Invalid command line option\n\n");
				return false;
			}
		}
	}

	// derive the other filenames
	char szbase[256]; strcpy(szbase, ops.szfile);
	char* ch = strrchr(szbase, '.');
	if (ch) *ch = 0;

	char* szext = (ch?ch+1:0);

	char szpath[256]; strcpy(szpath, ops.szfile);
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

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : Hello
// Prints the FEBio banner to a file
//

void Hello(FILE* fp)
{
	fprintf(fp,"===========================================================================\n");
	fprintf(fp,"         ________    _________   _________     __     _________            \n");
	fprintf(fp,"        |        |\\ |        |\\ |        |\\   |  |\\  /         \\\\          \n");
	fprintf(fp,"        |    ____|| |    ____|| |    __  ||   |__|| |    ___    ||         \n");
	fprintf(fp,"        |   |\\___\\| |   |\\___\\| |   |\\_| ||    \\_\\| |   //  \\   ||         \n");
	fprintf(fp,"        |   ||      |   ||      |   || | ||    __   |  ||    |  ||         \n");
	fprintf(fp,"        |   ||__    |   ||__    |   ||_| ||   |  |\\ |  ||    |  ||         \n");
	fprintf(fp,"        |       |\\  |       |\\  |         \\\\  |  || |  ||    |  ||         \n");
	fprintf(fp,"        |    ___||  |    ___||  |    ___   || |  || |  ||    |  ||         \n");
	fprintf(fp,"        |   |\\__\\|  |   |\\__\\|  |   |\\__|  || |  || |  ||    |  ||         \n");
	fprintf(fp,"        |   ||      |   ||      |   ||  |  || |  || |  ||    |  ||         \n");
	fprintf(fp,"        |   ||      |   ||___   |   ||__|  || |  || |   \\\\__/   ||         \n");
	fprintf(fp,"        |   ||      |        |\\ |          || |  || |           ||         \n");
	fprintf(fp,"        |___||      |________|| |__________|| |__||  \\_________//          \n");
	fprintf(fp,"                                                                           \n");
	fprintf(fp,"      F I N I T E   E L E M E N T S   F O R   B I O M E C H A N I C S      \n");
	fprintf(fp,"                                                                           \n");
	fprintf(fp,"                   --- v e r s i o n - %d . %d . %d ---                    \n", VERSION, SUBVERSION, SUBSUBVERSION);
	fprintf(fp,"                                                                           \n");
	fprintf(fp,"                                                                           \n");
	fprintf(fp,"  Musculoskeletal Research Laboratory                                      \n");
	fprintf(fp,"  University of Utah                                                       \n");
	fprintf(fp,"  http://mrl.sci.utah.edu                                                  \n");
	fprintf(fp,"                                                                           \n");
	fprintf(fp,"  copyright (c) 2006-2008 - All rights reserved                            \n");
	fprintf(fp,"                                                                           \n");
	fprintf(fp,"===========================================================================\n");
	fprintf(fp,"\n\n");
}

//-----------------------------------------------------------------------------

void init_framework(FEM& fem)
{
	FEBioCommand::SetFEM(&fem);
}
