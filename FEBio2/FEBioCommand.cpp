#include "stdafx.h"
#include <cstdlib>
#include "FEBioCommand.h"
#include "FEBioLib/FEBioModel.h"
#include "FEBioLib/version.h"
#include "FECore/FEException.h"
#include "FECore/FECore.h"
#include "NumCore/CompactMatrix.h"
#include "FECore/FEAnalysis.h"

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
REGISTER_COMMAND(FEBioCmd_Time   , "time"   , "print progress time statistics");

//-----------------------------------------------------------------------------


FEBioModel* FEBioCommand::m_pfem = 0;

FEBioCommand::FEBioCommand()
{
}

FEBioCommand::~FEBioCommand()
{
}

void FEBioCommand::SetFEM(FEBioModel* pfem)
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
	FEAnalysis* pstep = m_pfem->GetCurrentStep();
	bool bdebug = m_pfem->GetDebugFlag();
	if (nargs == 1) bdebug = !bdebug;
	else
	{
		if (strcmp(argv[1], "on") == 0) bdebug = true;
		else if (strcmp(argv[1], "off") == 0) bdebug = false;
		else { fprintf(stderr, "%s is not a valid option for debug.\n", argv[1]); return 0; }
	}
	m_pfem->SetDebugFlag(bdebug);

	printf("Debug mode is %s\n", (bdebug?"on":"off"));
	return 0;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Dtmin::run(int nargs, char **argv)
{
	assert(m_pfem);
	if (nargs == 2)
	{
		m_pfem->GetCurrentStep()->m_dtmin = atof(argv[1]);
		printf("Minumum time step size = %lg\n", m_pfem->GetCurrentStep()->m_dtmin);
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
	assert(false);
	return 1;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Print::run(int nargs, char **argv)
{
	assert(m_pfem);
	FEAnalysis* pstep = m_pfem->GetCurrentStep();

	if (nargs >= 2)
	{
		if (strcmp(argv[1], "time") == 0)
		{
			printf("Time : %lg\n", m_pfem->GetCurrentTime());
		}
		else
		{
			// assume it is a material parameter
			double* pd = m_pfem->FindParameter(ParamString(argv[1]));
			if (pd)
			{
				printf("%lg\n", *pd);
			}
			else
			{
				printf("The variable %s is not recognized\n", argv[1]);
			}
		}
	}
	else printf("Incorrect number of arguments for print command\n");

	return 0;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Restart::run(int nargs, char **argv)
{
	assert(m_pfem);
	FEAnalysis* pstep = m_pfem->GetCurrentStep();
	int ndump = pstep->GetDumpLevel();

	if (nargs == 1) ndump = (ndump == FE_DUMP_NEVER ? FE_DUMP_MAJOR_ITRS : FE_DUMP_NEVER);
	else
	{
		if      (strcmp(argv[1], "on" ) == 0) ndump = FE_DUMP_NEVER;
		else if (strcmp(argv[2], "off") == 0) ndump = FE_DUMP_MAJOR_ITRS;
		else 
		{
			fprintf(stderr, "%s is not a valid option for restart.\n", argv[1]);
			return 0;
		}
	}

	pstep->SetDumpLevel(ndump);

	return 0;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Version::run(int nargs, char **argv)
{
#ifdef _WIN64
	printf("\nFEBio version %d.%d.%d (x64)\n", VERSION, SUBVERSION, SUBSUBVERSION);
#else
	printf("\nFEBio version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
#endif
	printf("\nSVN revision: %d\n", SVNREVISION);
	printf("compiled on " __DATE__ "\n");
	printf("using FECore version %s\n\n", FECore::get_version_string());
	return 0;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Time::run(int nargs, char **argv)
{
	double sec = m_pfem->GetSolveTimer().peek();
	double sec0 = sec;

	int nhour, nmin, nsec;

	nhour = (int) (sec / 3600.0); sec -= nhour*3600;
	nmin  = (int) (sec /   60.0); sec -= nmin*60;
	nsec  = (int) (sec);
	printf("Elapsed time       :  %d:%02d:%02d\n", nhour, nmin, nsec);

	double endtime = m_pfem->GetCurrentStep()->m_tend;

	double pct = (m_pfem->GetCurrentTime() - m_pfem->GetCurrentStep()->m_dt) / endtime;
	if ((pct != 0) && (m_pfem->GetCurrentStep()->m_ntimesteps != 0))
	{
		double sec1 = sec0*(1.0/pct - 1.0);
		nhour = (int) (sec1 / 3600.0); sec1 -= nhour*3600;
		nmin  = (int) (sec1 /   60.0); sec1 -= nmin*60;
		nsec  = (int) (sec1);
		printf("Est. time remaining:  %d:%02d:%02d\n", nhour, nmin, nsec);
	}
	else
		printf("Est. time remaining:  (not available)\n");

	return 0;
}
