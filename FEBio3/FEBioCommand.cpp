#include "stdafx.h"
#include <cstdlib>
#include <FEBioLib/FEBioModel.h>
#include <FEBioLib/version.h>
#include <FECore/FEException.h>
#include <FECore/FESolver.h>
#include <FECore/CompactMatrix.h>
#include <FECore/FEAnalysis.h>
#include <NumCore/MatrixTools.h>
#include "FEBioCommand.h"
#include <iostream>
#include <fstream>

//-----------------------------------------------------------------------------
REGISTER_COMMAND(FEBioCmd_Cont         , "cont"   , "continues run");
REGISTER_COMMAND(FEBioCmd_Conv         , "conv"   , "force conversion of iteration");
REGISTER_COMMAND(FEBioCmd_Debug        , "debug"  , "toggle debug mode");
REGISTER_COMMAND(FEBioCmd_Fail         , "fail"   , "force iteratoin failer");
REGISTER_COMMAND(FEBioCmd_Help         , "help"   , "print available commands");
REGISTER_COMMAND(FEBioCmd_Plot         , "plot"   , "store current state to plot file");
REGISTER_COMMAND(FEBioCmd_Print        , "print"  , "print values of variables");
REGISTER_COMMAND(FEBioCmd_Quit         , "quit"   , "terminate the run and quit");
REGISTER_COMMAND(FEBioCmd_Version      , "version", "print version information");
REGISTER_COMMAND(FEBioCmd_Time         , "time"   , "print progress time statistics");
REGISTER_COMMAND(FEBioCmd_svg          , "svg"    , "write matrix sparsity pattern to svg file");
REGISTER_COMMAND(FEBioCmd_out          , "out"    , "write matrix and rhs file");
REGISTER_COMMAND(FEBioCmd_where        , "where"  , "current callback event");
REGISTER_COMMAND(FEBioCmd_break        , "break"  , "add a break point");
REGISTER_COMMAND(FEBioCmd_breaks       , "breaks" , "print list of break points");
REGISTER_COMMAND(FEBioCmd_clear_breaks , "clear"  , "clear one or all break points");

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
class ExitRequest : public std::runtime_error
{
public:
	ExitRequest() throw() : std::runtime_error("Early termination by user request") {}
};

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
			FEParamValue val = m_pfem->GetParameterValue(ParamString(argv[1]));
			if (val.isValid())
			{
				switch (val.type())
				{
				case FE_PARAM_DOUBLE: printf("%lg\n", val.value<double>()); break;
				default:
					printf("(cannot print value)\n");
				}
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

int FEBioCmd_Version::run(int nargs, char **argv)
{
#ifdef _DEBUG
	printf("\nFEBio version %d.%d.%d (DEBUG)\n", VERSION, SUBVERSION, SUBSUBVERSION);
#else
	printf("\nFEBio version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
#endif
	printf("compiled on " __DATE__ "\n");
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

	double endtime = m_pfem->GetEndTime();

	FETimeInfo& tp = m_pfem->GetTime();
	double pct = (tp.currentTime - tp.timeIncrement) / endtime;
	if (pct != 0)
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

//-----------------------------------------------------------------------------
int FEBioCmd_svg::run(int nargs, char **argv)
{
	FESolver* solver = m_pfem->GetCurrentStep()->GetFESolver();
	SparseMatrix* M = solver->GetStiffnessMatrix()->GetSparseMatrixPtr();
	std::vector<double> R = solver->GetLoadVector();
	CompactMatrix* A = dynamic_cast<CompactMatrix*>(M);
	if (A && m_pfem->GetFileTitle())
	{
		int rows = A->Rows();
		int cols = A->Columns();

		int i0 = 0, j0 = 0;
		int i1 = -1, j1 = -1;
		if (nargs == 3)
		{
			i1 = atoi(argv[1]);
			if (i1 < 0) { i0 = rows + i1; i1 = rows - 1; }
			j1 = atoi(argv[2]);
			if (j1 < 0) { j0 = cols + j1; j1 = cols - 1; }
		}

		const char* szfile = m_pfem->GetFileTitle();
		char buf[1024] = { 0 }, szsvg[1024] = { 0 };
		strcpy(buf, szfile);
		char* ch = strrchr(buf, '.');
		if (ch) *ch = 0;
		sprintf(szsvg, "%s.svg", buf);

		std::filebuf fb;
		fb.open(szsvg, std::ios::out);
		std::ostream out(&fb);
		NumCore::print_svg(A, out, i0, j0, i1, j1);
		fb.close();

		cout << "\nFile written " << szsvg << endl;
	}

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_out::run(int nargs, char **argv)
{
	FESolver* solver = m_pfem->GetCurrentStep()->GetFESolver();
	SparseMatrix* M = solver->GetStiffnessMatrix()->GetSparseMatrixPtr();
	std::vector<double> R = solver->GetLoadVector();
	CompactMatrix* A = dynamic_cast<CompactMatrix*>(M);
	if (A && m_pfem->GetFileTitle())
	{
		const char* szfile = m_pfem->GetFileTitle();
		char buf[1024] = { 0 }, szK[1024] = { 0 }, szR[1024] = { 0 };
		strcpy(buf, szfile);
		char* ch = strrchr(buf, '.');
		if (ch) *ch = 0;
		sprintf(szK, "%s.out", buf);
		sprintf(szR, "%s_rhs.out", buf);

		NumCore::write_hb(*A, szK);
		NumCore::write_vector(R, szR);

		cout << "\nFiles written: " << szK << ", " << szR << endl;
	}

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_where::run(int nargs, char **argv)
{
	unsigned int nevent = m_pfem->CurrentEvent();
	if (nevent == 0) cout << "Not inside callback event\n";
	else cout << "Callback event: ";

	switch (nevent)
	{
	case CB_INIT         : cout << "INIT"; break;
	case CB_STEP_ACTIVE  : cout << "STEP_ACTIVE"; break;
	case CB_MAJOR_ITERS  : cout << "MAJOR_ITERS"; break;
	case CB_MINOR_ITERS  : cout << "MINOR_ITERS"; break;
	case CB_SOLVED       : cout << "SOLVED"; break;
	case CB_UPDATE_TIME  : cout << "UPDATE_TIME"; break;
	case CB_AUGMENT      : cout << "AUGMENT"; break;
	case CB_STEP_SOLVED  : cout << "STEP_SOLVED"; break;
	case CB_MATRIX_REFORM: cout << "MATRIX_REFORM"; break;
	default:
		cout << "(unknown)";
	}
	cout << endl;

	return 0;
}

//-----------------------------------------------------------------------------

// in FEBio.cpp
void add_break_point(double t);
void add_cb_break_point(int nwhen);
void print_break_points();
void clear_break_points(int n);

int FEBioCmd_break::run(int nargs, char **argv)
{
	if (nargs != 2)
	{
		cout << "Invalid number of arguments";
	}

	const char* szbuf = argv[1];

	if      (_stricmp(szbuf, "ALWAYS"       ) == 0) add_cb_break_point(CB_ALWAYS);
	else if (_stricmp(szbuf, "INIT"         ) == 0) add_cb_break_point(CB_INIT);
	else if (_stricmp(szbuf, "STEP_ACTIVE"  ) == 0) add_cb_break_point(CB_STEP_ACTIVE);
	else if (_stricmp(szbuf, "MAJOR_ITERS"  ) == 0) add_cb_break_point(CB_MAJOR_ITERS);
	else if (_stricmp(szbuf, "MINOR_ITERS"  ) == 0) add_cb_break_point(CB_MINOR_ITERS);
	else if (_stricmp(szbuf, "SOLVED"       ) == 0) add_cb_break_point(CB_SOLVED);
	else if (_stricmp(szbuf, "UPDATE_TIME"  ) == 0) add_cb_break_point(CB_UPDATE_TIME);
	else if (_stricmp(szbuf, "AUGMENT"      ) == 0) add_cb_break_point(CB_AUGMENT);
	else if (_stricmp(szbuf, "STEP_SOLVED"  ) == 0) add_cb_break_point(CB_STEP_SOLVED);
	else if (_stricmp(szbuf, "MATRIX_REFORM") == 0) add_cb_break_point(CB_MATRIX_REFORM);
	else
	{
		double f = atof(szbuf);
		add_break_point(f);
	}

	return 0;
}

int FEBioCmd_breaks::run(int nargs, char **argv)
{
	print_break_points();
	return 0;
}

int FEBioCmd_clear_breaks::run(int nargs, char **argv)
{
	if (nargs == 1)
		clear_break_points(-1);
	else if (nargs == 2)
	{
		int bp = atoi(argv[1]);
		clear_break_points(bp - 1);
	}
	return 0;
}
