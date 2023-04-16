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
#include <cstdlib>
#include <FEBioLib/FEBioModel.h>
#include <FEBioLib/version.h>
#include <FECore/FEException.h>
#include <FECore/FESolver.h>
#include <FECore/CompactMatrix.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEGlobalMatrix.h>
#include "FEBioCommand.h"
#include "console.h"
#include <FEBioLib/cmdoptions.h>
#include <FEBioLib/febio.h>
#include <FEBioLib/plugin.h>
#include "FEBioApp.h"
#include "breakpoint.h"
#include <iostream>
#include <fstream>

//-----------------------------------------------------------------------------
// TODO: On Windows the GetCurrentTime macro gets in here via plugin.h. 
// I need to look into how to prevent this
#ifdef GetCurrentTime
#undef GetCurrentTime
#endif

//-----------------------------------------------------------------------------
#ifdef WIN32
#define szcmp    _stricmp
#else
#define szcmp    strcmp
#endif

//-----------------------------------------------------------------------------
REGISTER_COMMAND(FEBioCmd_break        , "break"  , "add a break point");
REGISTER_COMMAND(FEBioCmd_breaks       , "breaks" , "print list of break points");
REGISTER_COMMAND(FEBioCmd_clear_breaks , "clear"  , "clear one or all break points");
REGISTER_COMMAND(FEBioCmd_Config       , "config" , "(re-)load a FEBio configuration file");
REGISTER_COMMAND(FEBioCmd_Cont         , "cont"   , "continues the current model");
REGISTER_COMMAND(FEBioCmd_Conv         , "conv"   , "force conversion of iteration");
REGISTER_COMMAND(FEBioCmd_Debug        , "debug"  , "toggle debug mode");
REGISTER_COMMAND(FEBioCmd_Events       , "events" , "print list of events");
REGISTER_COMMAND(FEBioCmd_Fail         , "fail"   , "force iteratoin failer");
REGISTER_COMMAND(FEBioCmd_Help         , "help"   , "print available commands");
REGISTER_COMMAND(FEBioCmd_hist         , "hist"   , "lists history of commands");
REGISTER_COMMAND(FEBioCmd_LoadPlugin   , "import" , "load a plugin");
REGISTER_COMMAND(FEBioCmd_Plot         , "plot"   , "store current state to plot file");
REGISTER_COMMAND(FEBioCmd_out          , "out"    , "write matrix and rhs file");
REGISTER_COMMAND(FEBioCmd_Plugins      , "plugins", "list the plugins that are loaded");
REGISTER_COMMAND(FEBioCmd_Print        , "print"  , "print values of variables");
REGISTER_COMMAND(FEBioCmd_Quit         , "quit"   , "terminate the run and quit");
REGISTER_COMMAND(FEBioCmd_Restart      , "restart", "toggle restart mode");
REGISTER_COMMAND(FEBioCmd_Run          , "run"    , "run an FEBio input file");
REGISTER_COMMAND(FEBioCmd_set          , "set"    , "set value of some model and config parameters");
REGISTER_COMMAND(FEBioCmd_svg          , "svg"    , "write matrix sparsity pattern to svg file");
REGISTER_COMMAND(FEBioCmd_Time         , "time"   , "print progress time statistics");
REGISTER_COMMAND(FEBioCmd_UnLoadPlugin , "unload" , "unload a plugin");
REGISTER_COMMAND(FEBioCmd_Version      , "version", "print version information");
REGISTER_COMMAND(FEBioCmd_where        , "where"  , "current callback event");
REGISTER_COMMAND(FEBioCmd_list         , "list"   , "list factory classes");

int need_active_model()
{
	printf("No active model.\n");
	return 0;
}

int model_already_running()
{
	printf("A model is running. You must stop the active model before running this command.\n");
	return 0;
}

int invalid_nr_args()
{
	printf("Invalid number of arguments.\n");
	return 0;
}

int unknown_args()
{
	printf("Unrecognized arguments.\n");
	return 0;
}


//-----------------------------------------------------------------------------
FEBioCommand::FEBioCommand()
{
}

FEBioCommand::~FEBioCommand()
{
}

FEBioModel* FEBioCommand::GetFEM()
{
	return FEBioApp::GetInstance()->GetCurrentModel();
}

//-----------------------------------------------------------------------------
int FEBioCmd_Run::run(int nargs, char** argv)
{
	FEBioModel* fem = GetFEM();
	if (fem) return model_already_running();

	FEBioApp* febio = FEBioApp::GetInstance();

	if (febio->ParseCmdLine(nargs, argv) == false) return 0;

	// run FEBio on the ops
	febio->RunModel();

	// reset the title after computation.
	Console* pShell = Console::GetHandle();
	pShell->SetTitle("FEBio4");

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_Restart::run(int nargs, char** argv)
{
	FEBioModel* fem = GetFEM();
	if (fem == nullptr) return need_active_model();

	int dumpLevel = fem->GetDumpLevel();
	if (dumpLevel == FE_DUMP_NEVER) fem->SetDumpLevel(FE_DUMP_MAJOR_ITRS);
	else fem->SetDumpLevel(FE_DUMP_NEVER);

	dumpLevel = fem->GetDumpLevel();
	printf("Restart level set to: ");
	switch (dumpLevel)
	{
	case FE_DUMP_NEVER      : printf("NEVER (0)\n"); break;
	case FE_DUMP_MAJOR_ITRS : printf("MAJOR_ITRS (1)\n"); break;
	case FE_DUMP_STEP       : printf("STEP (2)\n"); break;
	case FE_DUMP_MUST_POINTS: printf("MUST POINTS (3)\n"); break;
	default:
		printf("(unknown value)\n");
		break;
	}

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_LoadPlugin::run(int nargs, char** argv)
{
	FEBioModel* fem = GetFEM();
	if (fem) return model_already_running();

	if (nargs < 2) fprintf(stderr, "missing file name\n");
	else febio::ImportPlugin(argv[1]);

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_UnLoadPlugin::run(int nargs, char* argv[])
{	
	FEBioModel* fem = GetFEM();
	if (fem) return model_already_running();

	FEBioPluginManager* PM = FEBioPluginManager::GetInstance();
	if (PM == 0) return -1;

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

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_Version::run(int nargs, char** argv)
{
	char* szver = febio::getVersionString();

#ifdef _DEBUG
	fprintf(stderr, "\nFEBio version %s (DEBUG)\n", szver);
#else
	fprintf(stderr, "\nFEBio version %s\n", szver);
#endif
	fprintf(stderr, "SDK Version %d.%d\n", FE_SDK_MAJOR_VERSION, FE_SDK_SUB_VERSION);
	fprintf(stderr, "compiled on " __DATE__ "\n\n");

	return 0;
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
		int l = (int)strlen(szn);
		printf("\t%s ", szn);
		while (l++ - 15 < 0) putchar('.');
		const char* szd = (*it)->GetDescription();
		printf(" : %s\n", szd);
	}

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_Events::run(int nargs, char** argv)
{
	printf("\n");
	printf("\tALWAYS        : break on any event.\n");
	printf("\tINIT          : break after model initialization.\n");
	printf("\tSTEP_ACTIVE   : break after step activation.\n");
	printf("\tMAJOR_ITERS   : break after major iteration converged.\n");
	printf("\tMINOR_ITERS   : break after minor iteration.\n");
	printf("\tSOLVED        : break after model is solved.\n");
	printf("\tUPDATE_TIME   : break before time is incremented.\n");
	printf("\tAUGMENT       : break before augmentation.\n");
	printf("\tSTEP_SOLVED   : break after step is solved.\n");
	printf("\tMATRIX_REFORM : break after global matrix is reformed.\n");
	printf("\n");
	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_Config::run(int nargs, char* argv[])
{
	FEBioModel* fem = GetFEM();
	if (fem) return model_already_running();

	FEBioApp* feApp = FEBioApp::GetInstance();
	febio::CMDOPTIONS& ops = feApp->CommandOptions();

	if (nargs == 1)
	{
		feApp->Configure(ops.szcnf);
	}
	else if (nargs == 2)
	{
		sprintf(ops.szcnf, "%s", argv[1]);
		feApp->Configure(ops.szcnf);
	}
	else return invalid_nr_args();

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_Plugins::run(int nargs, char** argv)
{
	FEBioPluginManager* PM = FEBioPluginManager::GetInstance();
	if (PM == 0) return 0;

	int NP = PM->Plugins();
	if (NP == 0)
	{
		fprintf(stdout, "no plugins loaded\n");
	}
	else
	{
		for (int i = 0; i < NP; ++i)
		{
			const FEBioPlugin& pl = PM->GetPlugin(i);
			fprintf(stdout, "%%%d: %s\n", i + 1, pl.GetName());
		}
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
	if (GetFEM()) throw ExitRequest();
	return 1;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Cont::run(int nargs, char** argv)
{
	if (GetFEM() == nullptr) return need_active_model();
	return 1;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Conv::run(int nargs, char **argv)
{
	if (GetFEM() == nullptr) return need_active_model();
	throw ForceConversion();
	return 1;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Debug::run(int nargs, char** argv)
{
	FEBioModel* fem = GetFEM();
	if (fem == nullptr)
	{
		printf("No active model.\n");
		return 0;
	}

	FEAnalysis* pstep = fem->GetCurrentStep();
	int ndebug = fem->GetDebugLevel();
	if (nargs == 1) ndebug = (ndebug? 0 : 1);
	else
	{
		if      (strcmp(argv[1], "on" ) == 0) ndebug = 1;
		else if (strcmp(argv[1], "off") == 0) ndebug = 0;
		else { fprintf(stderr, "%s is not a valid option for debug.\n", argv[1]); return 0; }
	}
	fem->SetDebugLevel(ndebug);

	printf("Debug mode is %s\n", (ndebug?"on":"off"));
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
	assert(false);
	return 1;
}

//-----------------------------------------------------------------------------

int FEBioCmd_Print::run(int nargs, char **argv)
{
	FEBioModel* fem = GetFEM();
	if (fem == nullptr) return need_active_model();

	FEAnalysis* pstep = fem->GetCurrentStep();

	if (nargs >= 2)
	{
		if (strcmp(argv[1], "time") == 0)
		{
			printf("Time : %lg\n", fem->GetCurrentTime());
		}
		else
		{
			// assume it is a material parameter
			FEParamValue val = fem->GetParameterValue(ParamString(argv[1]));
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
	else invalid_nr_args();

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_Time::run(int nargs, char **argv)
{
	FEBioModel* fem = GetFEM();
	if (fem == nullptr) return need_active_model();

	double sec = fem->GetSolveTimer().peek();
	double sec0 = sec;

	int nhour, nmin, nsec;

	nhour = (int) (sec / 3600.0); sec -= nhour*3600;
	nmin  = (int) (sec /   60.0); sec -= nmin*60;
	nsec  = (int) (sec);
	printf("Elapsed time       :  %d:%02d:%02d\n", nhour, nmin, nsec);

	double endtime = fem->GetEndTime();

	FETimeInfo& tp = fem->GetTime();
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
	FEBioModel* fem = GetFEM();
	if (fem == nullptr) return need_active_model();

	FESolver* solver = fem->GetCurrentStep()->GetFESolver();
	SparseMatrix* M = solver->GetStiffnessMatrix()->GetSparseMatrixPtr();
	std::vector<double> R = solver->GetLoadVector();
	CompactMatrix* A = dynamic_cast<CompactMatrix*>(M);
	if (A && (fem->GetFileTitle().empty() == false))
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

		const char* szfile = fem->GetFileTitle().c_str();
		char buf[1024] = { 0 }, szsvg[1024] = { 0 };
		strcpy(buf, szfile);
		char* ch = strrchr(buf, '.');
		if (ch) *ch = 0;
		sprintf(szsvg, "%s.svg", buf);

		std::filebuf fb;
		fb.open(szsvg, std::ios::out);
		std::ostream out(&fb);
		febio::print_svg(A, out, i0, j0, i1, j1);
		fb.close();

		cout << "\nFile written " << szsvg << endl;
	}

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_out::run(int nargs, char **argv)
{
	FEBioModel* fem = GetFEM();
	if (fem == nullptr) return need_active_model();

	int mode = 0; // binary
	if (nargs > 1)
	{
		if (strcmp(argv[1], "-txt") == 0) mode = 1; // text 
	}

	FESolver* solver = fem->GetCurrentStep()->GetFESolver();
	SparseMatrix* M = solver->GetStiffnessMatrix()->GetSparseMatrixPtr();
	std::vector<double> R = solver->GetLoadVector();
	CompactMatrix* A = dynamic_cast<CompactMatrix*>(M);
	if (A && (fem->GetFileTitle().empty() == false))
	{
		const char* szfile = fem->GetFileTitle().c_str();
		char buf[1024] = { 0 }, szK[1024] = { 0 }, szR[1024] = { 0 };
		strcpy(buf, szfile);
		char* ch = strrchr(buf, '.');
		if (ch) *ch = 0;
		sprintf(szK, "%s.out", buf);
		sprintf(szR, "%s_rhs.out", buf);

		febio::write_hb(*A, szK, mode);
		febio::write_vector(R, szR, mode);

		cout << "\nFiles written: " << szK << ", " << szR << endl;
	}
	else cout << "ERROR: Don't know how to write matrix format.\n";

	return 0;
}

//-----------------------------------------------------------------------------
int FEBioCmd_where::run(int nargs, char **argv)
{
	FEBioModel* fem = GetFEM();
	if (fem == nullptr) return need_active_model();

	unsigned int nevent = fem->CurrentEvent();
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
int FEBioCmd_break::run(int nargs, char **argv)
{
	if (nargs != 2) return invalid_nr_args();

	const char* szbuf = argv[1];

	add_break_point(szbuf);

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

int FEBioCmd_hist::run(int nargs, char **argv)
{
	Console* console = Console::GetHandle();
	vector<string> hist = console->GetHistory();
	int i = 1;
	for (string& s : hist)
	{
		cout << "%" << i++ << ": " << s << endl;
	}

	return 0;
}

const char* super_id_to_string(SUPER_CLASS_ID superID)
{
	const char* szclass = FECoreKernel::SuperClassString(superID);
	if (szclass == nullptr) szclass = "(unknown)";
	return szclass;
}

int FEBioCmd_list::run(int nargs, char** argv)
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();

	const char* szmod  = 0;
	const char* sztype = 0;
	const char* szpat = 0;
	int nmax = 0;
	int lpat = 0;

	// process argumets
	std::cin.clear();
	int nc = 1;
	while (nc < nargs)
	{
		if (strcmp(argv[nc], "-m") == 0)
		{
			if (nargs == 2)
			{
				int mods = fecore.Modules();
				for (int i = 0; i < mods; ++i)
				{
					const char* szmod = fecore.GetModuleName(i);
					printf("%d: %s\n", i+1, szmod);
				}
				printf("\n");
				return 0;
			}

			if (nargs <= nc + 1) return invalid_nr_args();
			szmod = argv[++nc];
		}
		else if (strcmp(argv[nc], "-c") == 0)
		{
			if (nargs <= nc + 1) return invalid_nr_args();
			sztype = argv[++nc];
		}
		else if (strcmp(argv[nc], "-n") == 0)
		{
			if (nargs <= nc + 1) return invalid_nr_args();
			nmax = atoi(argv[++nc]);
		}
		else if (strcmp(argv[nc], "-s") == 0)
		{
			if (nargs <= nc + 1) return invalid_nr_args();
			szpat = argv[++nc];
			lpat = (int)strlen(szpat);
			if (lpat == 0) szpat = 0;
		}
		else return unknown_args();
		++nc;
	}

	int facs = fecore.FactoryClasses();
	int n = 1, m = 0;
	for (int i = 0; i < facs; ++i)
	{
		const FECoreFactory* fac = fecore.GetFactoryClass(i);
		if (fac == nullptr) { printf("%3d: %s\n", n++, "(null)"); m++; }
		else
		{
			SUPER_CLASS_ID superID = fac->GetSuperClassID();
			const char* szclass = super_id_to_string(superID);

			int moduleId = fac->GetModuleID();
			const char* szmodule = fecore.GetModuleNameFromId(moduleId);
			if ((szmod == 0) || (szmodule && (strcmp(szmodule, szmod) == 0)))
			{
				if ((sztype == 0) || (szcmp(szclass, sztype) == 0))
				{
					const char* facType = fac->GetTypeStr();
					if ((szpat == 0) || (strstr(facType, szpat)))
					{
						if (szmodule) printf("%3d: %s.%s [%s]\n", n++, szmodule, facType, szclass);
						else printf("%3d: %s [%s]\n", n++, facType, szclass);
						m++;
					}
				}
			}
		}

		if ((nmax != 0) && (m >= nmax))
		{
			char ch = 0;
			do {
				printf("Continue (y or n)?");
				std::cin.get(ch);
			} while ((ch != 'y') && (ch != 'n'));
			m = 0;
			if (ch == 'n') break;
		}
	}
	printf("\n");

	return 0;
}

int FEBioCmd_set::run(int nargs, char** argv)
{
	FEBioModel* fem = GetFEM();
	if (nargs == 1)
	{
		printf("output_negative_jacobians = %d\n", (NegativeJacobian::m_boutput ? 1 : 0));
		if (fem)
		{
			printf("print_model_params        = %d\n", (fem->GetPrintParametersFlag() ? 1 : 0));
			printf("show_warnings_and_errors  = %d\n", (fem->ShowWarningsAndErrors() ? 1 : 0));
		}
		return 0;
	}

	if (nargs != 3)
	{
		printf("insufficient arguments.");
		return 0;
	}

	int n = atoi(argv[2]);

	if (strcmp(argv[1], "output_negative_jacobians") == 0)
	{
		NegativeJacobian::m_boutput = (n != 0);
		printf("output_negative_jacobians = %d", n);
	}
	else if (fem && strcmp(argv[1], "print_model_params") == 0)
	{
		fem->SetPrintParametersFlag(n != 0);
		printf("print_model_params = %d", n);
	}
	else if (fem && strcmp(argv[1], "show_warnings_and_errors") == 0)
	{
		fem->ShowWarningsAndErrors(n != 0);
		printf("show_warnings_and_errors = %d", n);
	}
	else
	{
		printf("don't know \"%s\"", argv[1]);
	}

	return 0;
}
