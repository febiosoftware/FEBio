// FESolver.cpp: implementation of the FESolver class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FESolver.h"
#include "fem.h"
#include "FEException.h"
#include <stdlib.h>
#include <ctype.h>
#include "FECore/FECore.h"

///////////////////////////////////////////////////////////////////////////////
// FESolver Construction/Destruction
//

FESolver::FESolver(FEM& fem) : m_fem(fem), m_log(fem.m_log)
{
	// default values
	m_Rtol = 1e10;
	m_Dtol = 0.001;
	m_Etol = 0.01;
	m_LStol = 0.9;
	m_LSmin = 0.01;
	m_LSiter = 5;
	m_maxups = 10;
	m_maxref = 15;
	m_cmax   = 1e5;

	m_niter = 0;

	// Stiffness matrix and linear solver are allocated in Init()
	m_pK = 0;
	m_pM = 0;
	m_psolver = 0;
}

FESolver::~FESolver()
{
	delete m_pK;		// clean up stiffnes matrix data
	delete m_pM;		// clean up mass matrix data
	delete m_psolver;	// clean up linear solver data
}


///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolver::CreateStiffness
//  Creates the global stiffness matrix
//

bool FESolver::CreateStiffness(bool breset)
{
	// clean up the solver
	m_psolver->Destroy();

	// clean up the stiffness matrix
	m_pK->Clear();

	// create the stiffness matrix
	if (m_pK->Create(m_fem, breset) == false) 
	{
		m_log.printf("FATAL ERROR: An error occured while building the stiffness matrix\n\n");
		return false;
	}

	// Do the preprocessing of the solver
	m_SolverTime.start();
	{
		if (!m_psolver->PreProcess(*m_pK)) throw FatalError();
	}
	m_SolverTime.stop();

	// done!
	return true;
}


///////////////////////////////////////////////////////////////////////////////
// FESolver::interrupt
//  This is the interruption handler which is called after a CTRL+C interruption
//

void FESolver::interrupt()
{
	char szcmd[256] = {0};

	// ask for user input
	printf("\nEnter command. Type help for more information.");

	const int MAX_VARS = 32;
	char*	argv[MAX_VARS];
	int   nargs;

	while (1)
	{
		// print the command prompt
		printf("\n>");

		// you must flush the input buffer before using gets
		fflush(stdin);

		// get the command
		fgets(szcmd, 255, stdin);

		// fgets does not remove '\n' so we'll do it ourselves
		char* ch = strrchr(szcmd, '\n');
		if (ch) *ch = 0;

		// parse the arguments
		nargs = 0;
		int n = 0;
		int l = strlen(szcmd);
		ch = szcmd;
		for (int i=0; i<=l; ++i, ++ch)
		{
			if (!isspace(*ch) && (*ch != 0))
			{
				if (n==0) argv[nargs++] = ch;
				++n;
			}
			else 
			{
				if (n!=0)
				{
					argv[nargs-1][n] = 0;
					n = 0;
				}
			}
		}

		if (strcmp(argv[0], "help") == 0)
		{
			printf("\nCommand overview:\n");
			printf("\tdebug <n> - set debug level to n, where n is 0, 1 or 2.\n");
			printf("\tdtmin - set the minimum time step size\n");
			printf("\tdump - toggle dump flag\n");
			printf("\tcont - continues run\n");
			printf("\tconv - forces to converge the current time step\n");
			printf("\tfail - stop the current iteration and retry (if auto time stepper enabled)\n");
			printf("\tplot - dump current state to plot database and continue\n");
			printf("\trestart - toggle restart flag\n");
			printf("\tquit - exits the application\n");
			printf("\tversion - prints version information\n");
		}
		else if (strcmp(argv[0], "dtmin") == 0)
		{
			if (nargs == 2)
			{
				m_fem.m_pStep->m_dtmin = atof(argv[1]);
				printf("Minumum time step size = %lg\n", m_fem.m_pStep->m_dtmin);
			}
			else printf("invalid number of arguments for dtmin\n");
		}
		else if (strcmp(argv[0], "restart") == 0)
		{
			if (nargs == 1) m_fem.m_pStep->m_bDump = 1;
			else
			{
				int n = atoi(argv[1]);
				if (n==0) m_fem.m_pStep->m_bDump = false;
				else if (n==1) m_fem.m_pStep->m_bDump = true;
				else fprintf(stderr, "unrecognized value %d\n", n);
			}
		}
		else if (strcmp(argv[0], "quit") == 0) throw ExitRequest();
		else if (strcmp(argv[0], "conv") == 0) throw ForceConversion();
		else if (strcmp(argv[0], "fail") == 0) throw IterationFailure();
		else if (strcmp(argv[0], "cont") == 0) return;
		else if (strcmp(argv[0], "plot") == 0)
		{
			m_fem.m_plot.Write(m_fem);
			return;
		}
		else if (strncmp(argv[0], "debug", 5) == 0)
		{
			int n = atoi(argv[1]);
			if ((n==0) || (n==1)) 
			{
				m_fem.SetDebugFlag(n==1);
				if (n==1) m_fem.m_pStep->SetPlotLevel(FE_PLOT_MINOR_ITRS);					
			}
			return;
		}
		else if (strcmp(argv[0], "dump") == 0)
		{
			m_fem.m_pStep->m_bDump = !m_fem.m_pStep->m_bDump;
			printf("Dump flag is %s\n", (m_fem.m_pStep->m_bDump?"on":"off"));
		}
		else if (strcmp(argv[0], "version") == 0)
		{
			printf("\nFEBio version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
			printf("compiled on " __DATE__ "\n");
			printf("using FECore version %s\n\n", FECore::get_version_string());
		}
		else if (strncmp(argv[0], "print", 5) == 0)
		{
			if (nargs == 2)
			{
				if (strcmp(argv[1], "nnz") == 0)
				{
					int nnz = m_pK->NonZeroes();
					printf("Nonzeroes in stiffness matrix: %d\n", nnz);
				}
				else if (strcmp(argv[1], "neq") == 0)
				{
					printf("Number of equations: %d\n", m_fem.m_neq);
				}
				else if (strcmp(argv[1], "time") == 0)
				{
					printf("Time : %lg\n", m_fem.m_ftime);
				}
				else
				{
					printf("The variable %s is not recognized\n", argv[1]);
				}
			}
			else printf("Incorrect number of arguments for print command\n");
		}
		else
		{
			printf("The command '%s' was not understood\nType help for an overview of valid commands.", szcmd);
		}
	}
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FESolver::Serialize(Archive& ar)
{
	if (ar.IsSaving())
	{
		ar << m_Dtol << m_Etol << m_Rtol << m_LSmin << m_LStol << m_LSiter;
		ar << m_maxups;
		ar << m_maxref;
		ar << m_cmax;
	}
	else
	{
		ar >> m_Dtol >> m_Etol >> m_Rtol >> m_LSmin >> m_LStol >> m_LSiter;
		ar >> m_maxups;
		ar >> m_maxref;
		ar >> m_cmax;
	}
}
