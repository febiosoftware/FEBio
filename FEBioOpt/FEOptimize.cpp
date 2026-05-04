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
#include "FEOptimize.h"
#include "FEOptimizeData.h"
#include "FECore/FECoreKernel.h"
#include "FECore/log.h"
#include "FECore/Timer.h"
#include "FEBioReport.h"

//-----------------------------------------------------------------------------
#define VERSION 2
#define SUB_VERSION	0

//-----------------------------------------------------------------------------
FEOptimize::FEOptimize(FEModel* pfem) : FECoreTask(pfem), m_opt(pfem)
{

}

//-----------------------------------------------------------------------------
bool FEOptimize::Init(const char* szfile)
{
	char szversion[32] = { 0 };
	snprintf(szversion, sizeof(szversion), "version %d.%d", VERSION, SUB_VERSION);
	feLog("P A R A M E T E R   O P T I M I Z A T I O N   M O D U L E\n%s\n\n", szversion);

	// read the data from the xml input file
	if (m_opt.Input(szfile) == false) return false;

	// do initialization
	bool ret = m_opt.Init();

	if (ret == false)
	{
		feLogErrorEx(m_opt.GetFEModel(), "Failed to initialize the optimization data.");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEOptimize::Run()
{
	Timer timer(false);
	timer.start();

	// solve the problem
	bool bret = m_opt.m_status = m_opt.Solve();

	timer.stop();
	double elapsedTime = timer.GetTime();
	char sztime[64];
	Timer::time_str(elapsedTime, sztime);
	feLog("\n\tTotal elapsed time : %s (%lg sec)\n\n", sztime, elapsedTime);

	if (bret)
		feLog("\n\n N O R M A L   T E R M I N A T I O N\n\n");
	else 
		feLog("\n\n E R R O R   T E R M I N A T I O N\n\n");

	if (m_opt.m_createReport)
		BuildReport();

	return bret;
}

void FEOptimize::BuildReport()
{
	FEBioReport report;

	report.SetTitle("Parameter Optimization Report");
	report.SetStatus(m_opt.m_status ? 1 : 0);
	report.SetOptionsFile(m_opt.m_filename);

	FEReportSection& section = report.AddSection("Optimization Results");

	section.AddValue("Total iterations", m_opt.m_niter);
	section.AddValue("Final objective value", m_opt.minObj);
	section.AddValue("Final regression coef", m_opt.minR2);

	std::vector<std::string> paramNames(m_opt.InputParameters());
	std::vector<double> paramValues(m_opt.InputParameters());

	// report the parameters for the minimal value
	for (int i = 0; i < m_opt.InputParameters(); ++i)
	{
		FEInputParameter& var = *m_opt.GetInputParameter(i);
		paramNames[i] = var.GetName();
		paramValues[i] = m_opt.amin[i];
	}

	FEReportTable& paramTable = section.AddTable();
	paramTable.AddColumn("Parameter", paramNames);
	paramTable.AddColumn("Optimal Value", paramValues);

	section.AddTableView(paramTable).SetCaption("Optimal parameter values.");

	FEDataFitObjective* obj = dynamic_cast<FEDataFitObjective*>(&m_opt.GetObjective());
	if (obj)
	{
		int n = obj->Measurements();
		std::vector<double> x(n), y(n), d(n);
		obj->GetXValues(x);
		obj->GetMeasurements(d);
		obj->EvaluateFunctions(y);

		std::string x_name = "x";
		std::string y_name = "y";
		std::string d_name = "data"; // "data" from options file

		const FEDataParameter* src = dynamic_cast<const FEDataParameter*>(obj->GetDataSource());
		if (src)
		{
			x_name = src->GetOrdinateName();
			y_name = src->GetParameterName();
		}

		FEReportTable& table = section.AddTable();
		table.AddColumn(x_name, x);
		table.AddColumn(y_name, y);
		table.AddColumn(d_name, d);

		FEReportChart& chart = section.AddChart(FEReportChart::Line);
		chart.AddDataSeries("fit")
			.AddData(FEReportChart::X, table.id, x_name)
			.AddData(FEReportChart::Y, table.id, y_name)
			.AddData(FEReportChart::Y, table.id, d_name);
		chart.SetCaption("Comparison of fitted values and measurements.");

		section.AddTableView(table).SetCaption("Table of fitted values and measurements.");
	}

	// create the report's file name from the options file name (change the extension to .febr)
	std::string filename = m_opt.m_filename;
	size_t lastdot = filename.find_last_of(".");
	if (lastdot != std::string::npos)
		filename = filename.substr(0, lastdot);

	if (filename.empty())
		filename = "report";

	filename += ".febr";
	
	report.Write(filename);
	feLog("Report written to %s\n", filename.c_str());
}