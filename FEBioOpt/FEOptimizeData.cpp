#include "stdafx.h"
#include "FEOptimizeData.h"
#include "FELMOptimizeMethod.h"
#include "FEOptimizeInput.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>

//=============================================================================

//-----------------------------------------------------------------------------
FEModelParameter::FEModelParameter(FEModel* fem) : FEInputParameter(fem)
{
	m_pd = 0;
}

//-----------------------------------------------------------------------------
bool FEModelParameter::SetParameter(const string& paramName)
{
	SetName(paramName);

	// find the variable
	FEModel& fem = *GetFEModel();
	double* pd = fem.FindParameter(paramName.c_str());
	if (pd == 0) return false;

	// store the pointer to the parameter
	m_pd = pd;

	return true;
}

//-----------------------------------------------------------------------------
double FEModelParameter::GetValue()
{
	if (m_pd) return *m_pd;
	return 0;
}

//-----------------------------------------------------------------------------
bool FEModelParameter::SetValue(double newValue)
{
	if (m_pd == 0) return false;
	(*m_pd) = newValue;
	return true;
}

//=============================================================================

//-------------------------------------------------------------------
bool fecb(FEModel* pmdl, unsigned int nwhen, void* pd)
{
	// get the optimizaton data
	FEObjectiveFunction& obj = *((FEObjectiveFunction*)pd);

	// get the FEM data
	FEModel& fem = *obj.GetFEM();

	// get the current time value
	double time = fem.m_ftime;

	// evaluate the current reaction force value
	double value = *(obj.m_pd);

	// add the data pair to the loadcurve
	FELoadCurve& lc = obj.ReactionLoad();
	lc.Add(time, value);

	return true;
}

FEObjectiveFunction::FEObjectiveFunction(FEModel* fem) : m_fem(fem)
{
	
}

bool FEObjectiveFunction::Init()
{
	// make sure we have a model
	if (m_fem == 0) return false;

	// find all the parameters
	m_pd = m_fem->FindParameter(m_szname);
	if (m_pd == 0) return false;

	// register callback
	m_fem->AddCallback(fecb, CB_INIT | CB_MAJOR_ITERS, (void*) this);

	return true;
}

void FEObjectiveFunction::Reset()
{
	FELoadCurve& lc = ReactionLoad();
	lc.Clear();
}

double FEObjectiveFunction::Evaluate(vector<double>& y)
{
	FELoadCurve& lc = GetLoadCurve(m_nlc);
	int ndata = lc.Points();
	y.resize(ndata);

	double chisq = 0.0;
	FELoadCurve& rlc = ReactionLoad();
	felog.printf("               CURRENT        REQUIRED      DIFFERENCE\n");
	for (int i = 0; i<ndata; ++i)
	{
		double xi = lc.LoadPoint(i).time;
		double y0i = lc.LoadPoint(i).value;

		y[i] = rlc.Value(xi);
		double dy = (y[i] - y0i);
		chisq += dy*dy;
		felog.printf("%5d: %15.10lg %15.10lg %15lg\n", i + 1, y[i], y0i, fabs(y[i] - y0i));
	}
	felog.printf("objective value: %lg\n", chisq);

	return chisq;
}

//=============================================================================

//-----------------------------------------------------------------------------
FEOptimizeData::FEOptimizeData(FEModel& fem) : m_fem(fem), m_obj(&fem)
{
	m_pSolver = 0;
	m_pTask = 0;
	m_niter = 0;
}

//-----------------------------------------------------------------------------
FEOptimizeData::~FEOptimizeData(void)
{
	delete m_pSolver;
}

//-----------------------------------------------------------------------------
bool FEOptimizeData::Init()
{
	// allocate default optimization solver if none specified in input file
	if (m_pSolver == 0) m_pSolver = new FELMOptimizeMethod;

	// allocate default solver if none specified in input file
	if (m_pTask == 0) m_pTask = fecore_new<FECoreTask>(FETASK_ID, "solve", &m_fem);

	// do the initialization of the task
	if (m_pTask->Init(0) == false) return false;

	// initialize the objective function
	if (m_obj.Init() == false) return false;

	// don't plot anything
	m_fem.GetCurrentStep()->SetPlotLevel(FE_PLOT_NEVER);

	return true;
}

//-----------------------------------------------------------------------------
bool FEOptimizeData::Solve()
{
	// make sure we have a task that will solve the FE model
	if (m_pTask == 0) return false;

	// go for it!
	int NVAR = (int) m_Var.size();
	vector<double> amin(NVAR, 0.0);
	vector<double> ymin;
	double minObj = 0.0;
	bool bret = m_pSolver->Solve(this, amin, ymin, &minObj);
	if (bret)
	{
		felog.SetMode(Logfile::LOG_FILE_AND_SCREEN);

		felog.printf("\nP A R A M E T E R   O P T I M I Z A T I O N   R E S U L T S\n\n");

		felog.printf("\tFunction values:\n\n");
		for (int i=0; i<(int) ymin.size(); ++i)
			felog.printf("\t\t%15lg\n", ymin[i]);

		felog.printf("\tTotal iterations ........ : %15d\n\n", m_niter);
		felog.printf("\tFinal objective value ... : %15lg\n\n", minObj);
		felog.printf("\tOptimal parameters:\n\n");
		// report the parameters for the minimal value
		for (int i = 0; i<NVAR; ++i)
		{
			FEInputParameter& var = *GetInputParameter(i);
			string name = var.GetName();
			felog.printf("\t\t%-15s = %.16lg\n", name.c_str(), amin[i]);
		}
	}

	return bret;
}

//-----------------------------------------------------------------------------
//! Read the data from the input file
//!
bool FEOptimizeData::Input(const char *szfile)
{
	FEOptimizeInput in;
	if (in.Input(szfile, this) == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
bool FEOptimizeData::RunTask()
{
	if (m_pTask == 0) return false;
	return m_pTask->Run();
}

//-----------------------------------------------------------------------------
//! solve the FE problem with a new set of parameters
bool FEOptimizeData::FESolve(const vector<double>& a)
{
	// increase iterator counter
	m_niter++;

	// reset objective function data
	FEObjectiveFunction& obj = GetObjective();
	obj.Reset();

	// set the input parameters
	int nvar = InputParameters();
	if (nvar != (int)a.size()) return false;
	for (int i = 0; i<nvar; ++i)
	{
		FEInputParameter& var = *GetInputParameter(i);
		var.SetValue(a[i]);
	}

	// report the new values
	felog.SetMode(Logfile::LOG_FILE_AND_SCREEN);
	felog.printf("\n----- Iteration: %d -----\n", m_niter);
	for (int i = 0; i<nvar; ++i)
	{
		FEInputParameter& var = *GetInputParameter(i);
		string name = var.GetName();
		felog.printf("%-15s = %lg\n", name.c_str(), var.GetValue());
	}

	// reset the FEM data
	FEModel& fem = GetFEM();
	fem.Reset();

	// solve the FE problem
	felog.SetMode(Logfile::LOG_NEVER);
	bool bret = RunTask();
	felog.SetMode(Logfile::LOG_FILE_AND_SCREEN);

	return bret;
}
