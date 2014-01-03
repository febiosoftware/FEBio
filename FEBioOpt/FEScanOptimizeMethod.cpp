#pragma once

#include "stdafx.h"
#include "FEScanOptimizeMethod.h"
#include "FECore/Logfile.h"
#include "FECore/FECoreKernel.h"
#include "FEBio2/console.h"

BEGIN_PARAMETER_LIST(FEScanOptimizeMethod, FEOptimizeMethod)
	ADD_PARAMETERV(m_inc, FE_PARAM_DOUBLEV, 32, "inc");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// declared in dllmain.cpp
extern FECoreKernel* pFEBio;

static Logfile& GetLogfile() { return pFEBio->GetLogfile(); }
static Console& GetConsole() { return pFEBio->GetConsole(); }

//-----------------------------------------------------------------------------
// FEScanOptimizeMethod
//-----------------------------------------------------------------------------

void fecb(FEModel* pfem, void* pd);

bool FEScanOptimizeMethod::Solve(FEOptimizeData* pOpt)
{
	m_pOpt = pOpt;
	FEOptimizeData& opt = *pOpt;

	// set the variables
	int ma = opt.Variables();
	vector<double> a(ma);
	for (int i=0; i<ma; ++i)
	{
		OPT_VARIABLE& var = opt.Variable(i);
		a[i] = var.m_min;
	}

	// set the FEM callback function
	FEModel& fem = opt.GetFEM();
	fem.AddCallback(fecb, CB_MAJOR_ITERS, &opt);

	// set the data
	OPT_OBJECTIVE& obj = opt.GetObjective();
	FELoadCurve& lc = opt.GetLoadCurve(obj.m_nlc);
	int ndata = lc.Points();
	vector<double> x(ndata), y0(ndata), y(ndata);
	for (int i=0; i<ndata; ++i) 
	{
		x[i] = lc.LoadPoint(i).time;
		y0[i] = lc.LoadPoint(i).value;
	}

	opt.m_niter = 0;
	Logfile& log = GetLogfile();

	bool bdone = false;
	double fmin = 0.0;
	vector<double> amin;
	do
	{
		// solve the problem
		if (FESolve(x, a, y) == false) return false;

		// calculate objective function
		double fobj = 0.0;
		for (int i=0; i<ndata; ++i)
		{
			double dy = y[i] - y0[i];
			fobj += dy*dy;
		}
		log.printf("Objective value: %lg\n", fobj);

		if ((fmin == 0.0) || (fobj < fmin))
		{
			fmin = fobj;
			amin = a;
		}

		// adjust indices
		for (int i=0; i<ma; ++i)
		{
			OPT_VARIABLE& vi = opt.Variable(i);
			if (a[i] >= vi.m_max)
			{
				if (i<ma-1)
				{
					a[i+1] += m_inc[i+1];
					a[i  ] = vi.m_min;
				}
				else bdone = true;
			}
			else 
			{
				a[i] += m_inc[i];
				break;
			}
		}
	}
	while (!bdone);

	log.printf("\n-------------------------------------------------------\n");
	for (int i=0; i<ma; ++i) 
	{
		OPT_VARIABLE& var = opt.Variable(i);
		log.printf("%-15s = %lg\n", var.m_szname, amin[i]);
	}
	log.printf("Objective value: %lg\n", fmin);

	return true;
}

//-----------------------------------------------------------------------------
bool FEScanOptimizeMethod::FESolve(vector<double> &x, vector<double> &a, vector<double> &y)
{
	// get the optimization data
	FEOptimizeData& opt = *m_pOpt;

	// increase iterator counter
	opt.m_niter++;

	// get the FEM data
	FEModel& fem = opt.GetFEM();

	// reset reaction force data
	FELoadCurve& lc = opt.ReactionLoad();
	lc.Clear();

	// set the material parameters
	int nvar = opt.Variables();
	for (int i=0; i<nvar; ++i)
	{
		OPT_VARIABLE& var = opt.Variable(i);
		*(var.m_pd) = a[i];
	}

	// reset the FEM data
	fem.Reset();

	Logfile& log = GetLogfile();
	log.SetMode(Logfile::FILE_AND_SCREEN);
	log.printf("\n----- Iteration: %d -----\n", opt.m_niter);
	for (int i=0; i<nvar; ++i) 
	{
		OPT_VARIABLE& var = opt.Variable(i);
		log.printf("%-15s = %lg\n", var.m_szname, a[i]);
	}

	// solve the FE problem
	log.SetMode(Logfile::NEVER);
	Console& pwnd = GetConsole();
	pwnd.Deactivate();

	bool bret = fem.Solve();

	log.SetMode(Logfile::FILE_AND_SCREEN);
	if (bret)
	{
		FELoadCurve& rlc = opt.ReactionLoad();
		int ndata = x.size();
		if (m_print_level == PRINT_VERBOSE) log.printf("               CURRENT        REQUIRED      DIFFERENCE\n");
		for (int i=0; i<ndata; ++i) 
		{
			y[i] = rlc.Value(x[i]);
//			if (m_print_level == PRINT_VERBOSE) clog.printf("%5d: %15.10lg %15.10lg %15lg\n", i+1, y[i], m_y0[i], fabs(y[i] - m_y0[i]));
		}
	}

	return bret;
}

