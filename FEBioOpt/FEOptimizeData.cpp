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
#include "FEOptimizeData.h"
#include "FELMOptimizeMethod.h"
#include "FEOptimizeInput.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
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

	// see if we found the parameter
	if (val.isValid() == false)
	{
		feLogError("Cannot find parameter %s", name.c_str());
		return false;
	}

    switch (val.type()) {
        case FE_PARAM_DOUBLE:
        {
            // make sure we have a valid data pointer
            double* pd = (double*) val.data_ptr();
            if (pd == 0)
            {
                feLogError("Invalid data pointer for parameter %s", name.c_str());
                return false;
            }
            
            // store the pointer to the parameter
            m_pd = pd;
        }
            break;
            
        case FE_PARAM_VEC2D:
        {
            // make sure we have a valid data pointer
            vec2d* vd = (vec2d*) val.data_ptr();
            if (vd == 0)
            {
                feLogError("Invalid data pointer for parameter %s", name.c_str());
                return false;
            }
            // store the pointer to the parameter
            m_pd = &vd->y();
        }
            break;
            
       default:
        {
            feLogError("Invalid parameter type for parameter %s", name.c_str());
            return false;
        }
            break;
    }

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
FEOptimizeData::FEOptimizeData(FEModel* fem) : m_fem(fem)
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
	if (m_pSolver == 0) m_pSolver = new FELMOptimizeMethod(GetFEModel());

	// allocate default solver if none specified in input file
	if (m_pTask == 0) m_pTask = fecore_new<FECoreTask>("solve", m_fem);

	// don't plot anything
	for (int i = 0; i < m_fem->Steps(); ++i)
	{
		m_fem->GetStep(i)->SetPlotLevel(FE_PLOT_NEVER);
	}

	// do the initialization of the task
	GetFEModel()->BlockLog();
	if (m_pTask->Init(0) == false) return false;
	GetFEModel()->UnBlockLog();

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
		feLog("\nP A R A M E T E R   O P T I M I Z A T I O N   R E S U L T S\n\n");

		feLog("\tFunction values:\n\n");
		for (int i=0; i<(int) ymin.size(); ++i)
			feLog("\t\t%15lg\n", ymin[i]);

		// evaluate final regression coefficient
		vector<double> y0;
		m_obj->GetMeasurements(y0);
		double minR2 = m_obj->RegressionCoefficient(y0, ymin);

		feLog("\tTotal iterations ........ : %15d\n\n", m_niter);
		feLog("\tFinal objective value ... : %15lg\n\n", minObj);
        feLog("\tFinal regression coef ... : %15lg\n\n", minR2);
		feLog("\tOptimal parameters:\n\n");
		// report the parameters for the minimal value
		for (int i = 0; i<NVAR; ++i)
		{
			FEInputParameter& var = *GetInputParameter(i);
			string name = var.GetName();
			feLog("\t\t%-15s = %.16lg\n", name.c_str(), amin[i]);
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
	feLog("\n----- Iteration: %d -----\n", m_niter);
	for (int i = 0; i<nvar; ++i)
	{
		FEInputParameter& var = *GetInputParameter(i);
		string name = var.GetName();
		feLog("%-15s = %lg\n", name.c_str(), var.GetValue());
	}

	// reset the FEM data
	FEModel& fem = *GetFEModel();
	fem.BlockLog();
	fem.Reset();
	fem.UnBlockLog();

	// solve the FE problem
	fem.BlockLog();
	bool bret = RunTask();
	fem.UnBlockLog();

	return bret;
}
