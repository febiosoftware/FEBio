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
#include "FEParameterSweep.h"
#include <FECore/XMLReader.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
#include <FEBioPlot/FEBioPlotFile.h>
#include <filesystem>

BEGIN_FECORE_CLASS(FESweepParam, FECoreClass)
	ADD_PARAMETER(m_paramName, "param_name");
	ADD_PARAMETER(m_min, "min");
	ADD_PARAMETER(m_max, "max");
	ADD_PARAMETER(m_step, "step");
END_FECORE_CLASS();

FESweepParam::FESweepParam(FEModel* fem)
{
	m_min = m_max = m_step = 0.0;
	m_pd = nullptr;
}

FESweepParam::FESweepParam(const FESweepParam& p)
{
	m_min = p.m_min;
	m_max = p.m_max;
	m_step = p.m_step;
	m_paramName = p.m_paramName;
	m_pd = p.m_pd;
}

void FESweepParam::operator = (const FESweepParam& p)
{
	m_min = p.m_min;
	m_max = p.m_max;
	m_step = p.m_step;
	m_paramName = p.m_paramName;
	m_pd = p.m_pd;
}

void FESweepParam::SetValue(double v)
{
	*m_pd = v;
}

BEGIN_FECORE_CLASS(FEParameterSweep, FECoreStudy)
	ADD_PROPERTY(m_params, "param")->SetDefaultType("parameter_sweep_param");
END_FECORE_CLASS();

FEParameterSweep::FEParameterSweep(FEModel* fem) : FECoreStudy(fem)
{
	m_niter = 0;
}

//! initialization
bool FEParameterSweep::Init()
{
	// disable the standard output from FEBio.
	DisableStandardOutput();

	// initialize the model
	if (GetFEModel()->Init() == false) return false;

	// check the parameters
	if (InitParams() == false) return false;

	// create the output plot file
	plotFileName = "parameter_sweep.xplt";
	std::string optionsFile = GetOptionsFile();
	if (!optionsFile.empty())
	{
		// replace the extension with .xplt
		plotFileName = std::filesystem::path(optionsFile).stem().string() + ".xplt";
	}
	xplt = new FEBioPlotFile(GetFEModel());
	xplt->Open(plotFileName.c_str());
	xplt->Write(0);

	return FECoreStudy::Init();
}

bool FEParameterSweep::InitParams()
{
	// Make sure we have something to do
	if (m_params.empty()) return false;

	// check the parameters
	FEModel& fem = *GetFEModel();
	for (size_t i = 0; i<m_params.size(); ++i)
	{
		FESweepParam& p = m_params[i];

		// find the variable
		string name = p.m_paramName;
		FEParamValue val = fem.GetParameterValue(ParamString(name.c_str()));

		// see if we found the parameter
		if (val.isValid() == false)
		{
			feLogError("Cannot find parameter %s", name.c_str());
			return false;
		}

		// see if it's the correct type
		if (val.type() != FE_PARAM_DOUBLE)
		{
			feLogError("Invalid parameter type for parameter %s", name.c_str());
			return false;
		}

		// make sure we have a valid data pointer
		double* pd = (double*)val.data_ptr();
		if (pd == 0)
		{
			feLogError("Invalid data pointer for parameter %s", name.c_str());
			return false;
		}

		// store the pointer to the parameter
		p.m_pd = pd;
	}

	return true;
}

//! Run the optimization module
bool FEParameterSweep::Execute()
{
	size_t ma = m_params.size();
	vector<double> a(ma);
	for (size_t i = 0; i<ma; ++i)
	{
		FESweepParam& pi = m_params[i];
		a[i] = pi.m_min;
	}

	// run the parameter sweep
	bool bdone = false;
	do
	{
		// solve the problem with the new input parameters
		if (FESolve(a) == false) return false;

		// update indices
		for (size_t i = 0; i<ma; ++i)
		{
			FESweepParam& pi = m_params[i];
			a[i] += pi.m_step;
			if (a[i] <= pi.m_max) break;
			else if (i<ma - 1) a[i] = pi.m_min;
			else { bdone = true; }
		}
	}
	while (!bdone);

	xplt->Close();

	return true;
}

bool FEParameterSweep::FESolve(const vector<double>& a)
{
	++m_niter;
	feLog("\n----- Iteration: %d -----\n", m_niter);

	// set the input parameters
	Run run;
	size_t nvar = m_params.size();
	assert(nvar == a.size());
	for (int i = 0; i<nvar; ++i)
	{
		FESweepParam& var = m_params[i];
		var.SetValue(a[i]);
		run.params.push_back(a[i]);

		string name = var.m_paramName;
		feLog("%-15s = %lg\n", name.c_str(), a[i]);
	}

	// write the current state to the plot file
	xplt->Write(m_niter);

	// reset the FEM data
	FEModel& fem = *GetFEModel();
	fem.BlockLog();

	// reset model
	fem.Reset();

	// solve the FE problem
	bool bret = fem.Solve();
	run.success = bret;
	m_runs.push_back(run);

	fem.UnBlockLog();

	return bret;
}

void FEParameterSweep::BuildReport(FEBioReport& report)
{
	FEReportSection& section = report.AddSection("Runs");
	FEReportTable& table = section.AddTable();

	int nparams = (int)m_params.size();
	std::vector<std::string> iterData(m_runs.size());
	// add the first column for the iteration number
	for (size_t j = 0; j < m_runs.size(); ++j)
		iterData[j] = std::string("#") + std::to_string(j + 1);
	table.AddColumn("Run", iterData);

	// add columns for each parameter
	std::vector<double> data(m_runs.size());
	for (int i = 0; i < nparams; ++i)
	{
		for (size_t j = 0; j < m_runs.size(); ++j)
			data[j] = m_runs[j].params[i];
		table.AddColumn(m_params[i].m_paramName, data);
	}
	std::vector<std::string> successData(m_runs.size());
	for (size_t j = 0; j < m_runs.size(); ++j)
		successData[j] = m_runs[j].success ? "Success" : "Failed";
	table.AddColumn("Status", successData);

	section.AddText("The following table shows the parameters and status for each run of the parameter sweep.");
	section.AddTableView(table)
		.SetCaption("Parameter sweep runs.");

	FEReportSection& outSec = report.AddSection("Output");
	outSec.AddText("The following files were generated during the parameter sweep:");
	outSec.AddFile(plotFileName, "Plot file containing the results of the parameter sweep");
}