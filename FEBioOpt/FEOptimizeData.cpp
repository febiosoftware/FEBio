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
bool FEModelParameter::Init()
{
	// find the variable
	FEModel& fem = *GetFEModel();
	string name = GetName();
	FEParamValue val = fem.GetParameterValue(ParamString(name.c_str()));
	if (val.isValid() == false) return false;
	if (val.type() != FE_PARAM_DOUBLE) return false;
	double* pd = (double*) val.data_ptr();
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

//-----------------------------------------------------------------------------
FEOptimizeData::FEOptimizeData(FEModel& fem) : m_fem(fem)
{
	m_pSolver = 0;
	m_pTask = 0;
	m_niter = 0;
	m_obj = 0;
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
	if (m_pTask == 0) m_pTask = fecore_new<FECoreTask>("solve", &m_fem);

	// do the initialization of the task
	if (m_pTask->Init(0) == false) return false;

	// initialize all input parameters
	for (int i=0; i<(int)m_Var.size(); ++i)
	{
		FEInputParameter* p = m_Var[i];
		if ((p==0) || (p->Init() == false)) return false;

		// set the initial value
		p->SetValue(p->InitValue());
	}

	// initialize the objective function
	if (m_obj == 0) return false;
	if (m_obj->Init() == false) return false;

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

/*		felog.printf("\tFunction values:\n\n");
		for (int i=0; i<(int) ymin.size(); ++i)
			felog.printf("\t\t%15lg\n", ymin[i]);
*/
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
