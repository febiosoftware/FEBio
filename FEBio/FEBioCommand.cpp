#include "StdAfx.h"
#include "FEBioCommand.h"
#include "fem.h"
#include "FEException.h"
#include "FECore/FECore.h"

//-----------------------------------------------------------------------------
REGISTER_COMMAND(FEBioCmd_Cont   , "cont"   , "continues run");
REGISTER_COMMAND(FEBioCmd_Conv   , "conv"   , "force conversion of iteration");
REGISTER_COMMAND(FEBioCmd_Debug  , "debug"  , "toggle debug mode");
REGISTER_COMMAND(FEBioCmd_Dtmin  , "dtmin"  , "set min time step size");
REGISTER_COMMAND(FEBioCmd_Fail   , "fail"   , "force iteratoin failer");
REGISTER_COMMAND(FEBioCmd_Help   , "help"   , "print available commands");
REGISTER_COMMAND(FEBioCmd_Plot   , "plot"   , "store current state to plot file");
REGISTER_COMMAND(FEBioCmd_Print  , "print"  , "print values of variables");
REGISTER_COMMAND(FEBioCmd_Quit   , "quit"   , "terminate the run and quit");
REGISTER_COMMAND(FEBioCmd_Restart, "restart", "toggles restart flag");
REGISTER_COMMAND(FEBioCmd_Version, "version", "print version information");

//-----------------------------------------------------------------------------


FEM* FEBioCommand::m_pfem = 0;

FEBioCommand::FEBioCommand()
{
}

FEBioCommand::~FEBioCommand()
{
}

void FEBioCommand::SetFEM(FEM* pfem)
{
	m_pfem = pfem;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Help::run(int nargs, char** argv)
{
	CommandManager* pCM = CommandManager::GetInstance();
	int N = pCM->Size();
	if (N == 0) return 0;

	printf("\nCommand overview:\n");

	CommandManager::CmdIterator it = pCM->First();
	for (int i=0; i<N; ++i, ++it)
	{
		const char* szn = (*it)->GetName();
		const char* szd = (*it)->GetDescription();
		printf("\t%s - %s\n", szn, szd);
	}

	return 0;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Quit::run(int nargs, char** argv)
{
	throw ExitRequest();
	return 1;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Cont::run(int nargs, char** argv)
{
	return 1;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Conv::run(int nargs, char **argv)
{
	throw ForceConversion();
	return 1;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Debug::run(int nargs, char** argv)
{
	assert(m_pfem);
	bool bdebug = m_pfem->GetDebugFlag();
	if (nargs == 1) bdebug = !bdebug;
	else
	{
		if (strcmp(argv[1], "on") == 0) bdebug = true;
		else if (strcmp(argv[1], "off") == 0) bdebug = false;
		else { fprintf(stderr, "%s is not a valid option for debug.\n", argv[1]); return 0; }
	}
	m_pfem->SetDebugFlag(bdebug);
	if (bdebug) m_pfem->m_pStep->SetPlotLevel(FE_PLOT_MINOR_ITRS); 
	else m_pfem->m_pStep->SetPlotLevel(FE_PLOT_MAJOR_ITRS);

	printf("Debug mode is %s\n", (bdebug?"on":"off"));
	return 0;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Dtmin::run(int nargs, char **argv)
{
	assert(m_pfem);
	if (nargs == 2)
	{
		m_pfem->m_pStep->m_dtmin = atof(argv[1]);
		printf("Minumum time step size = %lg\n", m_pfem->m_pStep->m_dtmin);
	}
	else printf("invalid number of arguments for dtmin\n");
	return 0;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Fail::run(int nargs, char **argv)
{
	throw IterationFailure();
}

//-----------------------------------------------------------------------------

int FEBioCmd_Plot::run(int nargs, char **argv)
{
	assert(m_pfem);
	m_pfem->m_plot.Write(*m_pfem);
	return 1;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Print::run(int nargs, char **argv)
{
	assert(m_pfem);
	FESolver* psolver = m_pfem->m_pStep->m_psolver;

	if (nargs == 2)
	{
		if (strcmp(argv[1], "nnz") == 0)
		{
			int nnz = psolver->m_pK->NonZeroes();
			printf("Nonzeroes in stiffness matrix: %d\n", nnz);
		}
		else if (strcmp(argv[1], "neq") == 0)
		{
			printf("Number of equations: %d\n", m_pfem->m_neq);
		}
		else if (strcmp(argv[1], "time") == 0)
		{
			printf("Time : %lg\n", m_pfem->m_ftime);
		}
		else
		{
			printf("The variable %s is not recognized\n", argv[1]);
		}
	}
	else printf("Incorrect number of arguments for print command\n");

	return 0;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Restart::run(int nargs, char **argv)
{
	assert(m_pfem);
	bool bdump = m_pfem->m_pStep->m_bDump;

	if (nargs == 1) bdump = !bdump;
	else
	{
		if (strcmp(argv[1], "on") == 0) bdump = true;
		else if (strcmp(argv[2], "off") == 0) bdump = false;
		else 
		{
			fprintf(stderr, "%s is not a valid option for restart.\n", argv[1]);
			return 0;
		}
	}

	m_pfem->m_pStep->m_bDump = bdump;
	printf("Restart flag is %s\n", (bdump?"on":"off"));

	return 0;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Version::run(int nargs, char **argv)
{
	printf("\nFEBio version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
	printf("compiled on " __DATE__ "\n");
	printf("using FECore version %s\n\n", FECore::get_version_string());
	return 0;
}
